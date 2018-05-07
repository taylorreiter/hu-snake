#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(Biostrings)
mh <- readAAStringSet(args[1])
uni <- readAAStringSet(args[2])
writeXStringSet(x = unique(c(uni, mh)), filepath = args[3], format = "fasta")
