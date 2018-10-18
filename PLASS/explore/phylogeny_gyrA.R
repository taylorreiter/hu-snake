# phylogeny of gyrA

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

info <- read.csv("inputs/hu_info.csv") # hu sample info
info <- filter(info, sample_origin == "SB1") # only keep SB1 samples

plass <- import_kegg(file = "outputs/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_SRR1976948.[0-9_]{2,12})", origin = "plass")

plass <-  filter(plass, bin %in% info$name) # filter to SB1 samples

# bin <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_[0-9]*)", origin = "bin")
# bin <- filter(bin %in% info$sample_origin)
# kegg_df <- rbind(crumb_df, bin_df) # combine two dfs
# kegg_df <- filter(kegg_df, bin %in% crumb_df$bin) # retain only relevant bins

# find marker genes ------------------------------------------------

# define crumb marker KO
gyrA <- 'K02469'
# gyrB <- 'K02470'
# recA <- 'K03553'

plass_gyrA <- filter(plass, geneID %in% gyrA) 
dim(plass_gyrA)
table(plass_gyrA$bin)
# hu-genome23 (Clostridiales bacterium 38_11) has 20 gyrA seqs, use these to build a quick phylogeny before making one for hu-genome36 (264) or hu-genome41 (476).
plass_gyrA <- filter(plass_gyrA, bin == "hu-genome23")
write.table(plass_gyrA[ , 1], "hu-genome23-gyrA-names.txt", row.names = F, col.names = F, quote = F)
# Bash ---------
# samtools faidx outputs/hu-genomes-plass-clean/sb1-plass.nbhd.nostop.99.fas
# while read inline 
# do
# samtools faidx outputs/hu-genomes-plass-clean/sb1-plass.nbhd.nostop.99.fas $inline >> hu-genome23-gyrA.fas
# done < hu-genome23-gyrA-names.txt


# new approach ------------------------------------------------------------

source("http://www.bioconductor.org/biocLite.R")
biocLite("msa")

library(seqinr)  
library(msa)
library(ape)

### Get clostridiales uniprot sequences into R. 
# 10 proteins found manually by searching the Uniprot database ( all firmicutes)
accNo <- c("AC=P94605 OR AC=P94604 OR AC=P0AES4 OR AC=Q45066 OR AC=P0C1U9 OR AC=O50628 OR AC=Q9X3Y3 OR AC=D8K235 OR AC=D8K235 OR AC=P35925")

# using the seqinr package
choosebank()  # show's available data bases
mybank <- choosebank(bank = "swissprot")
gyra_seq <- query("gyrASeq", accNo) # N.B. protein info NOT returned in the same order as requested
gyra_seqs <- getSequence(gyra_seq) # Extract sequences of interest
# check the names of the sequences
getName(gyra_seq)

write.fasta(sequences = gyra_seqs, 
            names = getName(gyra_seq),
            nbchar = 80, file.out = "gyra-firmicutes-seqs.faa") # put in fasta, which is needed for multialignment
#Bash -----------
# cat gyra-firmicutes-seqs.faa hu-genome23-gyrA.fas > gyra-hu23-firmicutes.faa
# -----
mySeqs <- readAAStringSet("gyra-hu23-firmicutes.faa")   # from package 
mySeqs <- mySeqs[-11]
mySeqs
#Biostrings
myAln <- msa(mySeqs) # Perform a multiple sequence alignment 

# Turn your alignment into a tree
myAln2 <- msaConvert(myAln, type="seqinr::alignment") # convert the alignment for seqinr
d <- dist.alignment(myAln2, "identity") # generate a distance matrix using seqinr (Fitch matrix)
myTree <- nj(d) # the nj() function allows neighbor-joining tree estimation in ape

# plot the tree
plot(myTree, main="Phylogenetic Tree of hu-genome 24 (Clostridiales) gyrA Sequences")

