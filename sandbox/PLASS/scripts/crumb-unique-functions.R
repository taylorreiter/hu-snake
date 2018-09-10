# are there any functions that appear in the crumb assemblies that do not appear at all in the query genomes?

library(dplyr)

# Read in all bin and crumb annotations -----------------------------------

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

crumb_df <- import_kegg(file = "outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_SRR1976948.[0-9_]{2,12})", origin = "crumb")
bin_df <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_[0-9]*)", origin = "bin")

bin_df <- bin_df %>% 
  filter(bin %in% crumb_df$bin)  %>% # retain only relevant bins 
  filter(is.na(bin) == F)
kegg_df <- rbind(crumb_df, bin_df) # combine two dfs



    
# apply to each bin -------------------------------------------------------
unique_crumb <- function(KEGG_df, Bin){
  library(dplyr)
  hu <- kegg_df %>% filter(bin == Bin)
  hu_crumb <- hu %>% filter(origin == "crumb")
  hu_bin <- hu %>% filter(origin == "bin")
  df<-hu_crumb[!(hu_crumb$geneID %in% hu_bin$geneID), ]
  return(df)
}
    
uniques <- lapply(unique(kegg_df$bin), function(x){unique_crumb(KEGG_df = kegg_df, Bin = x)})
names(uniques) <- unique(kegg_df$bin)
tmp <- do.call(rbind, uniques)

kegg_df[grep("K07540", kegg_df$geneID), ]

tmp <- kegg_df %>% filter(bin == 'hu-genome21') %>% filter(!is.na(geneID))
write.table(tmp$geneID, "~/Desktop/hu-genome21.txt", quote = F, row.names = F, col.names = F)
