library(clusterProfiler)
library(dplyr)

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

crumb_df <- import_kegg(file = snakemake@input[["crumb_kegg"]], 
                        gsub_regex = "(ass_[0-9]*)", origin = "crumb")

keggo <- crumb_df[ , 2] # grab kegg orthologs
keggo <- keggo[!is.na(keggo)] # remove NAs
kegg_enrich <- enrichKEGG(gene = keggo, organism = "ko") # enrichment analysis
plot5b <- dotplot(kegg_enrich) # plot enrichment

pdf(snakemake@output[["pdf"]], width = 7.5, height = 3.5)
plot5b
dev.off()

png(snakemake@output[["png"]], width = 7.5, height = 3.5, units = 'in', res = 300)
plot5b
dev.off()