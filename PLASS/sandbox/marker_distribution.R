library(dplyr)
library(rentrez)
library(Rsamtools)

setwd(dir = "~/github/hu-snake/sandbox/PLASS/")
# INPUT  ------------------------------------------------------------------

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

#info <- read.csv(snakemake@input[['info']]) # hu sample info
info <- read.csv("inputs/hu_info.csv")
info <- info %>% 
          filter(sample_origin == "SB1") 

# plass <- import_kegg(file = snakemake@input[['crumb_kegg']], gsub_regex = "(_SRR1976948.[0-9_]{2,12})", origin = "crumb")
plass <- import_kegg(file = "outputs/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_SRR1976948.[0-9_]{2,12})", origin = "crumb")

plass <- plass %>% 
              filter(bin %in% info$name) # filter to SB1 samples


# GRAB MARKER GENES -------------------------------------------------------

# define crumb marker KO
gyrA <- 'K02469'
# gyrB <- 'K02470'
# recA <- 'K03553'

plass_marker <- plass %>%
                  filter(geneID %in% c(gyrA)) #, gyrB, recA)) # find crumb marker genes

plass_marker$gene <- gsub("(; DNA gyrase subunit )([AB])( \\[EC:5.99.1.3\\])", "", plass_marker$name) # set gene name
# crumb_marker$gene <- gsub("; recombination protein RecA", "", crumb_marker$gene) # make a new column for gene symbol

write.table(plass_marker$locus_tag, "explore/prots/plass_gyra.txt", quote = F, row.names = F, col.names = F)
# WRITE OUT, DO SAMTOOLS FAIDX ON THE COMMAND LINE TO REDUCE FASTA SIZE.

## grab plass amino acid gyra sequences
## in bash:
# samtools faidx inputs/hu-genomes-plass100/sb1/all_sb1.nostop.c100.fas 
# while read inline
# do
#   samtools faidx inputs/hu-genomes-plass100/sb1/all_sb1.nostop.c100.fas $inline >> explore/prots/plass-gyra.fas
# done < explore/prots/plass_gyra.txt

## download ncbi protein archaea gyra sequences:
## files were created by navigating to NCBI and typing in "gyra" and limiting to either archaea or bacteria. The accession list was then downloaded. See here for more detailed instructions: https://taylorreiter.github.io/2018-09-05-ncbi-prot-download/

# while read inline 
# do
# i=$inline
# curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${i}&rettype=fasta&retmode=txt">>archaea-gyra.faa
# done < archaea-gyra.seq

# while read inline 
# do 
# i=$inline
# curl -s  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${i}&rettype=fasta&retmode=txt">>bacteria-gyra.faa
# done < bacteria-gyra.seq 

## Concatenate the sequences into a single prokaryotic gyrA seqs file
# cat archaea-gyra.faa bacteria-gyra.faa > prokaryotic-gyra.faa

## Make a blast database out of the sequences:
# makeblastdb -in prokaryotic-gyra.faa -dbtype prot -out prokaryotic-gyra-db/prokaryotic-gyra

## blast seqs
# blastp -db prokaryotic-gyra-db/prokaryotic-gyra -outfmt 6 -word_size 6 -gapopen 11 -gapextend 1 -window_size 40 -matrix BLOSUM62 -comp_based_stats 2 -max_target_seqs 1 -query plass-gyra.fas -out plass-gyra-blastp.tab
