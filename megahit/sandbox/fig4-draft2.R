library(tidyverse)
library(Biostrings)
library(Rsamtools)
library(rentrez)

setwd("~/github/hu-snake/")

info <- read.csv("inputs/hu_info.csv")
info <- info %>% 
          filter(sample_origin == "SB1")

blast <- read.csv("outputs/hu-crumbs-bin/assembly/marker-genes/blast_aa.csv")
import_kegg <- function(file, gsub_regex, origin){
  # file name: full path to user_ko_definition.txt
  # gsub_regex: regex to use on locus_tag column to create bin identifier. For crumbs, "(ass_[0-9]*)". For bins, "(_[0-9]*)"
  # origin: either crumb or bin
  df <- read.csv(file, sep = "\t", stringsAsFactors = FALSE, quote = "", na.strings = "", 
                 col.names = c("locus_tag", "geneID", "name", "score", "second", "second_score")) 
  df$bin <- gsub(gsub_regex, "", df$locus_tag) # add column named bin for bin of origin
  df$origin <- rep(origin, nrow(df)) # add row specifying identity
  df <- df[df$score > 100, ] # retain only kegg ids for which the score is > 100
  return(df)
}

crumb_df <- import_kegg(file = "outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt", gsub_regex = "(ass_[0-9]*)", origin = "crumb")
bin_df <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_[0-9]*)", origin = "bin")

# find marker genes ------------------------------------------------

# define crumb marker KO
gyrA <- 'K02469'
gyrB <- 'K02470'
recA <- 'K03553'

crumb_marker <- crumb_df %>%
                  filter(geneID %in% c(gyrA, gyrB, recA))  %>% # find crumb marker genes
                  filter(bin %in% info$name)

crumb_marker$gene <- gsub("(; DNA gyrase subunit )([AB])( \\[EC:5.99.1.3\\])", "", crumb_marker$name) # set gene name
crumb_marker$gene <- gsub("; recombination protein RecA", "", crumb_marker$gene) # make a new column for gene symbol

# find marker genes to perform pairwise alignments with
bin_marker <- bin_df %>%
                filter(geneID %in% c(gyrA, gyrB, recA)) %>%
                filter(bin %in% crumb_marker$bin & geneID %in% crumb_marker$geneID) # filter to marker genes that are found in crumb assemblies

# pairwise alignments between crumb and query marker genes ----------------

# get amino acid sequences from crumbs for pairwise assembly
indexFa("outputs/hu-crumbs-bin/assembly/prokka/all.faa") 
crumbfaa = FaFile("outputs/hu-crumbs-bin/assembly/prokka/all.faa")
crumbgraa = as(seqinfo(crumbfaa), "GRanges") 

get_aas <- function(seq_name, GRanges, FaFile){
  idxaa = pmatch(seq_name, names(GRanges)) 
  seqaa = getSeq(FaFile, GRanges[idxaa], as="AAStringSet") 
  return(seqaa)
}

crumbaas <- lapply(X = crumb_marker[ , 1], 
              FUN = function(x) get_aas(seq_name = x, GRanges = crumbgraa, FaFile = crumbfaa))


indexFa("outputs/hu-bins/prokka/all-bins.faa") 
binfaa = FaFile("outputs/hu-bins/prokka/all-bins.faa")
bingraa = as(seqinfo(binfaa), "GRanges") 

binaas <- lapply(X = bin_marker[ , 1], 
                   FUN = function(x) get_aas(seq_name = x, GRanges = bingraa, FaFile = binfaa))

# compare all records to get correct sets of pairwise alignments

df <- data.frame()
for(i in 1:nrow(crumb_marker)){
  for(j in 1:nrow(bin_marker)){
    if(crumb_marker$bin[i] == bin_marker$bin[j] & crumb_marker$geneID[i] == bin_marker$geneID[j]){
      a <- pairwiseAlignment(pattern = crumbaas[[i]], subject = binaas[[j]], type = "local")
      x <- nrow(df) + i
      df[x, 'query_crumb'] <- names(crumbaas[[i]])
      df[x, 'subject_bin'] <- names(binaas[[j]])
      df[x, 'gene'] <- crumb_marker$gene[i]
      df[x, 'nchar'] <- nchar(a)
      df[x, 'nedit'] <- nedit(a)
      df[x, 'nmatch'] <-nmatch(a)
      df[x, 'nmismatch'] <-nmismatch(a)
      df[x, 'pid'] <-pid(a)
      df[x, 'score'] <-score(a)
      df[x, 'type'] <-type(a)
    } else {
      next
    }
  }
}



df <- df %>% 
        filter(!is.na(query_crumb))

df <- df %>% 
        mutate(bin = gsub("ass_[0-9]{5,}", "", query_crumb)) %>%
        left_join(., info, by = c("bin" = "name")) %>%
        mutate(species = taxonomy) %>%
        select(query_crumb, nchar, pid, bin, gene, taxonomy, abbreviation, species) %>%
        filter(nchar > 10) 

# plot! -------------------------------------------------------------------

df <- df[order(df$pid, decreasing = T), ]
# df <-df %>% mutate(color = ifelse(origin == "blast", "red", "black"))

plt <- ggplot(df) + 
        geom_point(aes(x = pid, y = reorder(abbreviation, pid), col = gene, size = nchar)) +
        theme_bw()+
        scale_x_continuous(limits = c(87, 100), name = "Amino Acid Percent Identity") +
        labs(size = "amino acid length", y = "Query Bin", color = "marker gene")

plt

# dep: blast novel markers ---------------------------------------------------------

# get the leftovers to grab from blast
crumb_blast <- crumb_marker %>%
  filter(!(locus_tag %in% df$query_crumb))

aas <- lapply(X = crumb_blast[ , 1], 
              FUN = function(x) get_aas(seq_name = x, GRanges = crumbgraa, FaFile = crumbfaa))

blastf6 <- c('qseqid', 'sgi', 'sacc', 'pident', 'length', 'bitscore', 'evalue') # Fields for blast to return.
blast_aas <- data.frame('qseqid' = character(0),  'sgi' = character(0), 'sacc' = character(0), 'pident' = numeric(0), 
                        'length' = numeric(0), 'bitscore' = numeric(0), 'evalue' = numeric(0), 'query' = character(0))

for(i in 1:length(aas)){
  print(i)
  blastoutaa <- system2('.snakemake/conda/b63dce9d/bin/blastp',
                        c('-db', "'nr'", '-outfmt', sprintf('"6 %s"', paste(collapse=' ', blastf6)),
                          '-remote', '-word_size', "'6'", '-gapopen', "'11'", '-gapextend', "'1'", '-window_size',
                          "'40'", '-matrix', "'BLOSUM62'", '-comp_based_stats', '"2"', '-max_target_seqs', "'1'"),
                        input = as.character(aas[[i]][1]), stdout = TRUE)
  blastdfaa <- `names<-`(read.table(quote="", sep='\t', textConnection(blastoutaa)), blastf6)
  blastdfaa$query <- rep(names(aas[[i]]), nrow(blastdfaa))
  blast_aas <- rbind(blast_aas, blastdfaa)
}


block2 <- function(i){
  pid_species <- entrez_summary(db="protein", id = i)
  Sys.sleep(2) # sleep; don't overwhelm NCBI server
  return(pid_species)
}

format_blast <- function(blast_df, rentrez_function, gsub_regex, info_df, crumb_df){
  blast_df <- blast_df[order(blast_df$query, blast_df$bitscore), ] # retain only best blast hit
  blast_df <- blast_df[!duplicated(blast_df$query, fromLast=TRUE), ] # Select the last row by id
  id_names <- lapply(unique(blast_df$sacc), rentrez_function) # grab info for each blast id
  id_species <- sapply(id_names, "[[", 26)   # extract species names from retrieved data
  id_df <- as.data.frame(unique(as.character(blast_df$sacc))) # make a dataframe of nuid and species
  id_df$species <- id_species # set species name
  blast_df <- merge(blast_df, id_df, by.x = "sacc", by.y = "unique(as.character(blast_df$sacc))") # merge dataframe with nuid species
  blast_df$name <- gsub(gsub_regex, "", blast_df$query) # parse blast queries to bin name
  blast_df <- left_join(blast_df, info_df, by = "name")
  blast_df <- left_join(blast_df, crumb_df, by = c("query" = "locus_tag"))
  return(blast_df)
}

blast_aas <- format_blast(blast_df = blast_aas, rentrez_function = block2, gsub_regex = "(ass_)([0-9]{5})", 
                          info_df = info, crumb_df = crumb_blast)


tmp <- blast_aas # save results to protect overwrite of obj
