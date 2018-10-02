# 1. Find gyrA, gyrB, and recA sequences in the crumb assemblies by parsing the KEGG GhostKoala output
# 2. write gyrA, gyrB, and recA nucleotide sequences to a fasta file
# 3. write gyrA, gyrB, and recA amino acid sequences to a fasta file
# 4. BLAST nucleotide sequences against NCBI 
# 5. BLAST amino acid sequences against NCBI
# 6. Report in a table. 

library(dplyr)
library(rentrez)
library(Rsamtools)

# read in data ------------------------------------------------------------
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
info <- read.csv(snakemake@input[['info']]) # hu sample info
info <- info %>% 
  filter(sample_origin == "SB1") # only keep SB1 samples

crumb_df <- import_kegg(file = snakemake@input[['crumb_kegg']], gsub_regex = "(ass_[0-9]*)", origin = "crumb")
crumb_df <- crumb_df %>% 
  filter(bin %in% info$name) # filter to SB1 samples

bin_df <- import_kegg(file =snakemake@input[['bin_kegg']], gsub_regex = "(_[0-9]*)", origin = "bin")

kegg_df <- rbind(crumb_df, bin_df) # combine two dfs
kegg_df <- kegg_df %>% 
              filter(bin %in% crumb_df$bin) # retain only relevant bins
# find marker genes ------------------------------------------------

# define crumb marker KO
gyrA <- 'K02469'
gyrB <- 'K02470'
recA <- 'K03553'

crumb_marker <- kegg_df %>% 
  filter(origin == "crumb") %>%
  filter(geneID %in% c(gyrA, gyrB, recA)) # find crumb marker genes

crumb_marker$gene <- gsub("(; DNA gyrase subunit )([AB])( \\[EC:5.99.1.3\\])", "", crumb_marker$name) # set gene name
crumb_marker$gene <- gsub("; recombination protein RecA", "", crumb_marker$gene) # make a new column for gene symbol

# blast nucleotides ----------------------------------------------------

indexFa(snakemake@input[['ffn']])   # create an index of file
# indexFa("outputs/hu-crumbs-bin/assembly/prokka/all.ffn")
ffn = FaFile(snakemake@input[['ffn']])
# ffn = FaFile("outputs/hu-crumbs-bin/assembly/prokka/all.ffn")
gr = as(seqinfo(ffn), "GRanges")

get_nucs <- function(seq_name, GRanges, FaFile){
  idx = pmatch(seq_name, names(GRanges)) 
  seq = getSeq(FaFile, GRanges[idx])
  return(seq)
}

nucs <- lapply(X = crumb_marker[ , 1], FUN = function(x) get_nucs(seq_name = x, GRanges = gr, FaFile = ffn))

# https://support.bioconductor.org/p/86440/
blastf6 <- c('qseqid', 'sgi', 'sacc', 'pident', 'length', 'bitscore', 'evalue') # Fields for blast to return.

blast_nuc <- data.frame('qseqid' = character(0),  'sgi' = character(0), 'sacc' = character(0), 'pident' = numeric(0), 
                        'length' = numeric(0), 'bitscore' = numeric(0), 'evalue' = numeric(0), 'query' = character(0)) #  initiate empty data frame for for loop results.

for(i in 1:length(nucs)){
  if(file.exists(snakemake@output[['blast_nuc']])){
    next #  if file exists do not run for loop.
  }
  print(i)
  Sys.sleep(30) # sleep; don't overwhelm NCBI server
  blastout <- system2(snakemake@input[['blastn']],
                      c('-db', "'nt'", '-outfmt', sprintf('"6 %s"', paste(collapse=' ', blastf6)),
                        '-remote', '-word_size', "'11'", '-gapopen', "'5'", '-gapextend', "'2'", '-reward',
                        "'2'", '-penalty', "'-3'", '-max_target_seqs', "'10'"),
                      input = as.character(nucs[[i]][1]), stdout = TRUE) #  run blast via the command line but save STD out to object to blastout.
  blastdf <- `names<-`(read.table(quote="", sep='\t', textConnection(blastout)), blastf6) # parse blastout as a table and assign names to it
  blastdf$query <- rep(names(nucs[[i]]), nrow(blastdf)) # add query name to blastdf
  blast_nuc <- rbind(blast_nuc, blastdf)
}

# if file doesn't exist, write to doc to avoid costly NCBI queries.
if(!file.exists(snakemake@output[['blast_nuc']])) {
  write.csv(blast_nuc, snakemake@output[['blast_nuc']], quote = F, row.names = F)
} else {
  blast_nuc <- read.csv(snakemake@output[['blast_nuc']])
}

block <- function(i){
  nuid_species <- entrez_summary(db="nucleotide", id = i)
  Sys.sleep(2) # sleep; don't overwhelm NCBI server
  return(nuid_species)
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

blast_nuc <- format_blast(blast_df = blast_nuc, rentrez_function = block, gsub_regex = "(ass_)([0-9]{5})", 
                          info_df = info, crumb_df = crumb_marker)

write.csv(blast_nuc, snakemake@output[['blast_nuc_full']], quote = F, row.names = F)

# blast amino acids -------------------------------------------------------------

indexFa(snakemake@input[['faa']]) 
faa = FaFile(snakemake@input[['faa']])
graa = as(seqinfo(faa), "GRanges") 

get_aas <- function(seq_name, GRanges, FaFile){
  idxaa = pmatch(seq_name, names(GRanges)) 
  seqaa = getSeq(FaFile, GRanges[idxaa], as="AAStringSet") 
  return(seqaa)
}

aas <- lapply(X = crumb_marker[ , 1], 
              FUN = function(x) get_aas(seq_name = x, GRanges = graa, FaFile = faa))

blast_aas <- data.frame('qseqid' = character(0),  'sgi' = character(0), 'sacc' = character(0), 'pident' = numeric(0), 
                        'length' = numeric(0), 'bitscore' = numeric(0), 'evalue' = numeric(0), 'query' = character(0))

for(i in 1:length(aas)){
  if(file.exists(snakemake@output[['blast_aas']])){
    next
  }
  print(i)
  Sys.sleep(30) # sleep; don't overwhelm NCBI server
  blastoutaa <- system2(snakemake@input[['blastp']],
                        c('-db', "'nr'", '-outfmt', sprintf('"6 %s"', paste(collapse=' ', blastf6)),
                          '-remote', '-word_size', "'6'", '-gapopen', "'11'", '-gapextend', "'1'", '-window_size',
                          "'40'", '-matrix', "'BLOSUM62'", '-comp_based_stats', '"2"', '-max_target_seqs', "'10'"),
                        input = as.character(aas[[i]][1]), stdout = TRUE)
  
  blastdfaa <- `names<-`(read.table(quote="", sep='\t', textConnection(blastoutaa)), blastf6)
  blastdfaa$query <- rep(names(aas[[i]]), nrow(blastdfaa))
  blast_aas <- rbind(blast_aas, blastdfaa)
}

if(!file.exists(snakemake@output[['blast_aas']])) {
  write.csv(blast_aas, snakemake@output[['blast_aas']], quote = F, row.names = F)
} else {
  blast_aas <- read.csv(snakemake@output[['blast_aas']])
}

block2 <- function(i){
  pid_species <- entrez_summary(db="protein", id = i)
  Sys.sleep(2) # sleep; don't overwhelm NCBI server
  return(pid_species)
}

blast_aas <- format_blast(blast_df = blast_aas, rentrez_function = block2, gsub_regex = "(ass_)([0-9]{5})", 
                          info_df = info, crumb_df = crumb_marker)


write.csv(blast_aas, snakemake@output[['blast_aas_full']], quote = F, row.names = F)

# non-snakemake code ------------------------------------------------------

# 1. Find gyrA, gyrB, and recA sequences in the crumb assemblies by parsing the KEGG GhostKoala output
# 2. write gyrA, gyrB, and recA nucleotide sequences to a fasta file
# 3. write gyrA, gyrB, and recA amino acid sequences to a fasta file
# 4. BLAST nucleotide sequences against NCBI 
# 5. BLAST amino acid sequences against NCBI
# 6. Report in a table. 

# library(dplyr)
# library(rentrez)
# library(Rsamtools)
# 
# # read in data
# import_kegg <- function(file, gsub_regex, origin){
#   # file name: full path to user_ko_definition.txt
#   # gsub_regex: regex to use on locus_tag column to create bin identifier. For crumbs, "(ass_[0-9]*)". For bins, "(_[0-9]*)"
#   # origin: either crumb or bin
#   df <- read.csv(file, sep = "\t", stringsAsFactors = FALSE, quote = "", na.strings = "",
#                  col.names = c("locus_tag", "geneID", "name", "score", "second", "second_score"))
#   df$bin <- gsub(gsub_regex, "", df$locus_tag) # add column named bin for bin of origin
#   df$origin <- rep(origin, nrow(df)) # add row specifying identity
#   df <- df[df$score > 100, ] # retain only kegg ids for which the score is > 100
#   return(df)
# }
# info <- read.csv("inputs/hu_info.csv") # hu sample info
# info <- info %>%
#   filter(sample_origin == "SB1") # only keep SB1 samples
# 
# crumb_df <- import_kegg(file = "outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt", gsub_regex = "(ass_[0-9]*)", origin = "crumb")
# crumb_df <- crumb_df %>%
#   filter(bin %in% info$name) # filter to SB1 samples
# 
# bin_df <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_[0-9]*)", origin = "bin")
# 
# kegg_df <- rbind(crumb_df, bin_df) # combine two dfs
# kegg_df <- kegg_df %>%
#   filter(bin %in% crumb_df$bin) # retain only relevant bins
# # find marker genes
# 
# # define crumb marker KO
# gyrA <- 'K02469'
# gyrB <- 'K02470'
# recA <- 'K03553'
# 
# crumb_marker <- kegg_df %>%
#   filter(origin == "crumb") %>%
#   filter(geneID %in% c(gyrA, gyrB, recA)) # find crumb marker genes
# 
# crumb_marker$gene <- gsub("(; DNA gyrase subunit )([AB])( \\[EC:5.99.1.3\\])", "", crumb_marker$name) # set gene name
# crumb_marker$gene <- gsub("; recombination protein RecA", "", crumb_marker$gene) # make a new column for gene symbol
# 
# indexFa("outputs/hu-crumbs-bin/assembly/prokka/all.ffn")   # create an index of file
# ffn = FaFile("outputs/hu-crumbs-bin/assembly/prokka/all.ffn")
# gr = as(seqinfo(ffn), "GRanges")
# 
# get_nucs <- function(seq_name, GRanges, FaFile){
#   idx = pmatch(seq_name, names(GRanges))
#   seq = getSeq(FaFile, GRanges[idx])
#   return(seq)
# }
# 
# nucs <- lapply(X = crumb_marker[ , 1], FUN = function(x) get_nucs(seq_name = x, GRanges = gr, FaFile = ffn))
# 
# # https://support.bioconductor.org/p/86440/
# blastf6 <- c('qseqid', 'sgi', 'sacc', 'pident', 'length', 'bitscore', 'evalue') # Fields for blast to return.
# 
# blast_nuc <- data.frame('qseqid' = character(0),  'sgi' = character(0), 'sacc' = character(0), 'pident' = numeric(0),
#                         'length' = numeric(0), 'bitscore' = numeric(0), 'evalue' = numeric(0), 'query' = character(0)) #  initiate empty data frame for for loop results.
# 
# for(i in 1:length(nucs)){
#   if(file.exists('outputs/hu-crumbs-bin/assembly/marker-genes/blast_nuc.csv')){
#     next #  if file exists do not run for loop.
#   }
#   print(i)
#   Sys.sleep(2) # sleep; don't overwhelm NCBI server
#   blastout <- system2('.snakemake/conda/b63dce9d/bin/blastn',
#                       c('-db', "'nt'", '-outfmt', sprintf('"6 %s"', paste(collapse=' ', blastf6)),
#                         '-remote', '-word_size', "'11'", '-gapopen', "'5'", '-gapextend', "'2'", '-reward',
#                         "'2'", '-penalty', "'-3'", '-max_target_seqs', "'10'"),
#                       input = as.character(nucs[[i]][1]), stdout = TRUE) #  run blast via the command line but save STD out to object to blastout.
#   blastdf <- `names<-`(read.table(quote="", sep='\t', textConnection(blastout)), blastf6) # parse blastout as a table and assign names to it
#   blastdf$query <- rep(names(nucs[[i]]), nrow(blastdf)) # add query name to blastdf
#   blast_nuc <- rbind(blast_nuc, blastdf)
# }
# 
# # if file doesn't exist, write to doc to avoid costly NCBI queries.
# if(!file.exists('outputs/hu-crumbs-bin/assembly/marker-genes/blast_nuc.csv')) {
#   write.csv(blast_nuc, 'outputs/hu-crumbs-bin/assembly/marker-genes/blast_nuc.csv', quote = F, row.names = F)
# } else {
#   blast_nuc <- read.csv('outputs/hu-crumbs-bin/assembly/marker-genes/blast_nuc.csv')
# }
# 
# block <- function(i){
#   nuid_species <- entrez_summary(db="nucleotide", id = i)
#   Sys.sleep(2) # sleep; don't overwhelm NCBI server
#   return(nuid_species)
# }
# 
# format_blast <- function(blast_df, rentrez_function, gsub_regex, info_df, crumb_df){
#   blast_df <- blast_df[order(blast_df$query, blast_df$bitscore), ] # retain only best blast hit
#   blast_df <- blast_df[!duplicated(blast_df$query, fromLast=TRUE), ] # Select the last row by id
#   id_names <- lapply(unique(blast_df$sacc), rentrez_function) # grab info for each blast id
#   id_species <- sapply(id_names, "[[", 26)   # extract species names from retrieved data
#   id_df <- as.data.frame(unique(as.character(blast_df$sacc))) # make a dataframe of nuid and species
#   id_df$species <- id_species # set species name
#   blast_df <- merge(blast_df, id_df, by.x = "sacc", by.y = "unique(as.character(blast_df$sacc))") # merge dataframe with nuid species
#   blast_df$name <- gsub(gsub_regex, "", blast_df$query) # parse blast queries to bin name
#   blast_df <- left_join(blast_df, info_df, by = "name")
#   blast_df <- left_join(blast_df, crumb_df, by = c("query" = "locus_tag"))
#   return(blast_df)
# }
# 
# blast_nuc <- format_blast(blast_df = blast_nuc, rentrez_function = block, gsub_regex = "(ass_)([0-9]{5})",
#                           info_df = info, crumb_df = crumb_marker)
# 
# write.csv(blast_nuc, 'outputs/hu-crumbs-bin/assembly/marker-genes/blast_nuc_full.csv', quote = F, row.names = F)
# 
# # blast amino acids
# 
# indexFa("outputs/hu-crumbs-bin/assembly/prokka/all.faa")
# faa = FaFile("outputs/hu-crumbs-bin/assembly/prokka/all.faa")
# graa = as(seqinfo(faa), "GRanges")
# 
# get_aas <- function(seq_name, GRanges, FaFile){
#   idxaa = pmatch(seq_name, names(GRanges))
#   seqaa = getSeq(FaFile, GRanges[idxaa], as="AAStringSet")
#   return(seqaa)
# }
# 
# aas <- lapply(X = crumb_marker[ , 1],
#               FUN = function(x) get_aas(seq_name = x, GRanges = graa, FaFile = faa))
# 
# blast_aas <- data.frame('qseqid' = character(0),  'sgi' = character(0), 'sacc' = character(0), 'pident' = numeric(0),
#                         'length' = numeric(0), 'bitscore' = numeric(0), 'evalue' = numeric(0), 'query' = character(0))
# 
# for(i in 1:length(aas)){
#   if(file.exists("outputs/hu-crumbs-bin/assembly/marker-genes/blast_aa.csv")){
#     next
#   }
#   print(i)
#   Sys.sleep(10)
#   blastoutaa <- system2('.snakemake/conda/b63dce9d/bin/blastp',
#                         c('-db', "'nr'", '-outfmt', sprintf('"6 %s"', paste(collapse=' ', blastf6)),
#                           '-remote', '-word_size', "'6'", '-gapopen', "'11'", '-gapextend', "'1'", '-window_size',
#                           "'40'", '-matrix', "'BLOSUM62'", '-comp_based_stats', '"2"', '-max_target_seqs', "'10'"),
#                         input = as.character(aas[[i]][1]), stdout = TRUE)
# 
#   blastdfaa <- `names<-`(read.table(quote="", sep='\t', textConnection(blastoutaa)), blastf6)
#   blastdfaa$query <- rep(names(aas[[i]]), nrow(blastdfaa))
#   blast_aas <- rbind(blast_aas, blastdfaa)
# }
# 
# if(!file.exists('outputs/hu-crumbs-bin/assembly/marker-genes/blast_aa.csv')) {
#   write.csv(blast_aas, 'outputs/hu-crumbs-bin/assembly/marker-genes/blast_aa.csv', quote = F, row.names = F)
# } else {
#   blast_aas <- read.csv('outputs/hu-crumbs-bin/assembly/marker-genes/blast_aa.csv')
# }
# 
# block2 <- function(i){
#   pid_species <- entrez_summary(db="protein", id = i)
#   Sys.sleep(2) # sleep; don't overwhelm NCBI server
#   return(pid_species)
# }
# 
# blast_aas <- format_blast(blast_df = blast_aas, rentrez_function = block2, gsub_regex = "(ass_)([0-9]{5})",
#                           info_df = info, crumb_df = crumb_marker)
# 
# 
# write.csv(blast_aas, 'outputs/hu-crumbs-bin/assembly/marker-genes/blast_aas_full.csv', quote = F, row.names = F)
