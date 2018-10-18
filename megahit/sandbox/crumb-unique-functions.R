# are there any functions that appear in the crumb assemblies that do not appear at all in the query genomes?

setwd("/Users/taylorreiter/github/hu-snake")
library(tidyverse)

# Read in all bin and crumb annotations -----------------------------------

# crumbs
    crumb_full <- list.files("outputs/hu-crumbs_bin/assembly/GhostKOALA/", "full.txt$", full.names = T)
    crumb_kegg <- list.files("outputs/hu-crumbs_bin/assembly/GhostKOALA/", "full.txt$")
    crumb_kegg <- gsub(".ko-ann-full.txt", "", crumb_kegg)
    crumb_df <- data.frame(locus_tag = NA, geneID = NA, name = NA, score = NA, 
                           second = NA, second_score = NA, bin = NA)
    for(i in 1:(length(crumb_full))){
      kegg <- read.csv(crumb_full[i], sep = "\t", stringsAsFactors = FALSE, na.strings = "",
                       col.names = c("locus_tag", "geneID", "name", "score", "second", "second_score", "bin")) 
      crumb_df <- rbind(crumb_df, kegg)
    }
    crumb_df$origin <- rep("crumb", nrow(crumb_df)) # add row specifying identity

# bins
    bin_full <- list.files("sandbox/bin_plus_crumb_assembly/bin", "full.txt$", full.names = T)
    bin_full <- bin_full[2:31]
    bin_kegg <- list.files("sandbox/bin_plus_crumb_assembly/bin/", "full.txt$")
    bin_kegg <- bin_kegg[2:31]
    bin_kegg <- gsub(".ko-ann-full.txt", "", bin_kegg)
    bin_df <- data.frame(locus_tag = NA, geneID = NA, name = NA, score = NA, 
                         second = NA, second_score = NA, bin = NA)
    for(i in 1:(length(bin_full))){
      kegg <- read.csv(bin_full[i], sep = "\t", stringsAsFactors = FALSE, na.strings = "",
                       col.names = c("locus_tag", "geneID", "name", "score", "second", "second_score", "bin")) 
      bin_df <- rbind(bin_df, kegg)
    }
    bin_df$origin <- rep("bin", nrow(bin_df))

# bind two dataframes together
    kegg_df <- rbind(crumb_df, bin_df)

# retain only kegg ids for which the score is > 100
    kegg_df <- kegg_df[kegg_df$score > 100, ]

# read in brite data so it's easier to understand function
    brite <- read.delim("sandbox/ko00001_parse.txt", header= F, sep = "\t")
    brite <- separate(brite, V4, into = c("id", "description"), sep = " ", 
                      remove= T, extra = "merge", fill = "left")

    kegg_df <- left_join(x = kegg_df, y = brite, by = c("geneID" = "id"))
    
# apply to each bin -------------------------------------------------------
unique.crumb <- function(KEGG_df, Bin){
  library(dplyr)
  hu <- kegg_df %>% filter(bin == Bin)
  hu_crumb <- hu %>% filter(origin == "crumb")
  hu_bin <- hu %>% filter(origin == "bin")
  df<-hu_crumb[!(hu_crumb$geneID %in% hu_bin$geneID), ]
  return(df)
}
    
uniques <- lapply(unique(kegg_df$bin)[2:31], function(x){unique.crumb(KEGG_df = kegg_df, Bin = x)})
names(uniques) <- unique(kegg_df$bin)[2:31]
tmp <- do.call(rbind, uniques)
tmp$edited <- gsub("(\\s\\[*)([A-Z]*)(:ko\\d*\\])", "", tmp$V3)
tmp$edited <- gsub("(\\d*)", "", tmp$edited)

ggplot(tmp, aes(x = edited, fill = bin)) + geom_bar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_flip()

# apply to all bins -------------------------------------------------------

crumbs <- kegg_df %>% filter(origin == "crumbs")
bins <- kegg_df %>% filter(origin == "bin")
!(crumbs$geneID %in% bins$geneID)
