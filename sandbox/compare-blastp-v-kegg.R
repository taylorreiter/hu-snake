#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Compare kegg id assignments to blastp
# args[1]: parsed blastxml file
# args[2]: kegg ghost koala file
# args[3]: output file

# Note xml format preserved the prokka assignment next to the amino acid sequence identifier assigned by prokka. Therefore can compare kegg, prokka, and blast

blast <- read.delim(args[1], header = TRUE)
kegg <- read.delim(args[2], header = FALSE)

# split blast first column on space to get contig id by itself
blast <-  tidyr::separate(data = blast, col = Query_Name, into = c("contig", "prokka"),
                   sep = " ", remove = TRUE, extra = "merge")

# remove NAs
blast <- blast[!is.na(blast$Bitscore), ]

# Grab highest bitscore match 
blast <- blast[order(blast$contig, blast$Bitscore), ]

# Select the last row by id
blast_bitscore <- blast[!duplicated(blast$contig, fromLast=TRUE), ]

merged2 <- merge(blast_bitscore, kegg, by.x = "contig", by.y = "V1")
write.table(merged2, file = args[3], row.names = F, quote = F, sep = "\t")


print(paste0("parsed ", args[1], " and ", args[2], ", output is in ", args[3]))
