library(rentrez)
library(dplyr)
library(ape)
library(Biostrings)
library(Rsamtools)
library(stringr)

# READ FILES --------------------------------------------------------------
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

info <- read.csv("inputs/hu_info.csv")
plass <- import_kegg(file = "outputs99/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_SRR1976948.[0-9_]{2,12})", origin = "crumb")
# FORMAT BLAST OUTPUT -----------------------------------------------------
blast <- read.delim("explore/prots/plass-archaea-gyra-blastp.tab", sep = "\t", header = F, 
                    col.names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                  'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'))
blast$bin <- gsub("(_SRR1976948.[0-9_]{2,12})", "", blast$qseqid)

block2 <- function(i){
  pid_species <- entrez_summary(db="protein", id = i)
  Sys.sleep(2) # sleep; don't overwhelm NCBI server
  return(pid_species)
}

format_blast <- function(blast_df, rentrez_function, gsub_regex, info_df, crumb_df){
  blast_df <- blast_df[order(blast_df$qseqid, blast_df$bitscore), ] # retain only best blast hit
  blast_df <- blast_df[!duplicated(blast_df$qseqid, fromLast=TRUE), ] # Select the last row by id
  id_names <- lapply(unique(blast_df$sseqid), rentrez_function) # grab info for each blast id
  id_species <- sapply(id_names, "[[", 26)   # extract species names from retrieved data
  id_df <- as.data.frame(unique(as.character(blast_df$sseqid))) # make a dataframe of nuid and species
  id_df$species <- id_species # set species name
  blast_df <- merge(blast_df, id_df, by.x = "sseqid", by.y = "unique(as.character(blast_df$sseqid))") # merge dataframe with nuid species
  blast_df$name <- gsub(gsub_regex, "", blast_df$qseqid) # parse blast queries to bin name
  blast_df <- left_join(blast_df, info_df, by = "name")
  blast_df <- left_join(blast_df, crumb_df, by = c("qseqid" = "locus_tag"))
  return(blast_df)
}

blast <- format_blast(blast_df = blast, rentrez_function = block2, gsub_regex = "(_SRR1976948.[0-9_]{2,12})", info_df = info, crumb_df = plass)


# WRANGLE BLAST -----------------------------------------------------------

# For each bin, find the highest percent match. 
# Grab the accession number of that match, and download it. 
# pairwise align all sequences to the best match
# record percent identity, and plot it as a bean plot

tops <- blast %>% group_by(name.x) %>% top_n(1, pident) # find highest match
tops <- tops[!duplicated(tops$name.x, fromLast=TRUE), ] # reduce to one

# grab sequences from local fasta that blast database was made from
indexFa("explore/prots/archaea-gyra.faa") 
faa = FaFile("explore/prots/archaea-gyra.faa")
graa = as(seqinfo(faa), "GRanges") 

get_aas <- function(seq_name, GRanges, FaFile){
  idxaa = pmatch(seq_name, names(GRanges)) 
  seqaa = getSeq(FaFile, GRanges[idxaa], as="AAStringSet") 
  return(seqaa)
}

blastaas <- sapply(X = as.character(tops$sseqid), 
              FUN = function(x) get_aas(seq_name = x, GRanges = graa, FaFile = faa))
blastaas <- AAStringSet(sapply(blastaas, `[[`, 1)) 
names(blastaas) <- paste0(as.character(tops$sseqid), "_", tops$bin.x)

# writeXStringSet(blastaas, file ="explore/prots/best-archaea-gyra.fasta", format = "fasta", append = T)

## In bash:
# makeblastdb -in best-archaea-gyra.fasta -dbtype prot
# blastp -db best-archaea-gyra.fasta -outfmt 6 -word_size 6 -gapopen 11 -gapextend 1 -window_size 40 -matrix BLOSUM62 -comp_based_stats 2 -max_target_seqs 1 -query plass-archaea-gyra.fas -out best-archaea-gyra-blastp.tab

# PAIRWISE ALIGN ----------------------------------------------------------
plassaas <- readAAStringSet("explore/prots/plass-archaea-gyra.fas")
# blastaas <- readAAStringSet("explore/prots/best-archaea-gyra.fasta")

# compare all records to get correct sets of pairwise alignments
df <- data.frame()
pwa <- list()
for(i in 1:length(plassaas)){ # loop through all plass aas
  print(paste0("The amino acid sequence is plassaas ", i))
  print(paste0("The faa to pairwise align is ", names(plassaas)[[i]]))
  for(bin in unique(tops$name.x)){ # loop through bin names
    print(paste0(bin))
    if(!(str_detect(string = names(plassaas)[[i]], pattern = bin))){
      next
    } else {
      for(j in 1:length(blastaas)){ # loop through best blast matches
        print(paste0("check if pwa should be with best aa match to ", names(blastaas)[j]))
        if(str_detect(string = names(blastaas[j]), pattern = bin)){
          print(paste0("Reference ", names(blastaas[j]), "and query ", names(plassaas[i]), "selected!"))
        #if(names(blastaas[j]) == bin){
            a <- pairwiseAlignment(pattern = plassaas[[i]], subject = blastaas[j], gapOpening = 11, 
                                   gapExtension = 1, substitutionMatrix = "BLOSUM62", type = "local") 
            pwa[[i]] <- a 
            x <- i
            df[x, 'query_plass'] <- names(plassaas[i])
            df[x, 'best_blast'] <- names(blastaas[j])
            # df[x, 'gene'] <- crumb_marker$gene[i]
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
    } 
  }

#df <- df %>% 
#        left_join(., info, by = c("best_blast" = "name")) %>%
#        mutate(species = taxonomy) %>%
#        select(query_plass, nchar, pid, bin, gene, taxonomy, abbreviation, species) %>%
#        filter(nchar > 0) 
df <- mutate(df, bin = gsub("(_SRR1976948.[0-9_]{2,12})", "", query_plass))
df <- df %>% 
        left_join(., info, by = c("bin" = "name"))

# PLOT --------------------------------------------------------------------
library(beanplot)

par(mar=c(5, 13, 1, 1) +.1) #bottom, left, top, right
beanplot(df$pid ~ df$bin, ll = .04, side = "second", method = "stack", names = unique(as.character(df$abbreviation)), xlab = "Amino Acid Percent Identity by Pairwise Local Alignment", horizontal = T, show.names = T, las = 1)
