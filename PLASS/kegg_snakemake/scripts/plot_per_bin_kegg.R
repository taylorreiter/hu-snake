# this script visualizes the per-bin difference in functional content between the PLASS assembled crumbs and the bins, and analyzes the amount of unique content in the PLASS assembled crumbs.

setwd("/Users/taylorreiter/github/hu-snake/sandbox/PLASS/")

library(clusterProfiler)
library(dplyr)
library(ggplot2)

# READ & FORMAT DATA ------------------------------------------------------------

# read in metadata
info <- read.csv("inputs/hu_info.csv")

import_kegg <- function(file, gsub_regex, origin){
  # file name: full path to user_ko_definition.txt
  # gsub_regex: regex to use on locus_tag column to create bin identifier. For crumbs, "". For bins, "(_[0-9]*)"
  # origin: either crumb or bin
  df <- read.csv(file, sep = "\t", stringsAsFactors = FALSE, quote = "", na.strings = "", 
                 col.names = c("locus_tag", "geneID", "name", "score", "second", "second_score")) 
  df$bin <- gsub(gsub_regex, "", df$locus_tag) # add column named bin for bin of origin
  df$origin <- rep(origin, nrow(df)) # add row specifying identity
  df <- df[df$score > 100, ] # retain only kegg ids for which the score is > 100
  return(df)
}

crumb_df <- import_kegg(file = "outputs/GhostKOALA/user_ko_definition.txt", 
                        gsub_regex = "(_SRR1976948.[0-9_]{2,12})", origin = "crumb")

bin_df <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", 
                      gsub_regex = "(_[0-9]*)", origin = "bin")
bin_df <- bin_df %>%
            filter(bin %in% crumb_df$bin)

# combine two dfs
kegg_df <- rbind(crumb_df, bin_df)

# Get mappings of KO id to BRITE hierarchy description
keggo <- kegg_df[ , 2]
keggo <- keggo[!is.na(keggo)]

# set p value to 1 so everything will be returned;
# i.e. this is not actually an enrichment analysis
kegg_enrich <- enrichKEGG(gene = keggo, organism = "ko", pvalueCutoff = 1) 
dotplot(kegg_enrich)
# grab data.frame of results
df <- kegg_enrich@result

# PLOT ------------------------------------------------------------------

# Use kegg df to map original kegg occurrences to categories
mapping <- data.frame(Description = NA, geneID = NA)
for(i in 1:nrow(df)){
  tmp <- strsplit(x = df$geneID[i], split = "/")
  tmp2 <- rep(df$Description[i], times = length(tmp))
  tmp_df <- data.frame(Description = tmp2, geneID = tmp)
  colnames(tmp_df) <- c("Description", "geneID")
  mapping <- rbind(mapping, tmp_df) 
}

# Transfer mappings to original kegg data
transfer_mappings <- function(kegg_df, mappings, df_label){
  # kegg_df: df of tsv downloaded from Ghost Koala indicating locus_tag kegg id annotations
  # mappings: mappings of kegg ID to BRITE hierarchy (see above)
  # df_label: a string that will get added to the output df indicating which sample is in the df
  # i.e. "kegg19"
  library(dplyr)
  kegg_df <- right_join(kegg_df, mappings)
  kegg_df <- kegg_df[!is.na(kegg_df$locus_tag), ]
  kegg_df <- kegg_df[!is.na(kegg_df$geneID), ]
  kegg_df$sample <- rep(df_label, nrow(kegg_df))
  return(kegg_df)
}

# Transfer mappings
map <- transfer_mappings(kegg_df = kegg_df, mappings = mapping, df_label = "all")

# tally the description; use to set order of the plot
sum_df <- map %>%
            group_by(Description) %>%
            tally() 

# set order from overall sums
sum_df$Description <- factor(sum_df$Description, 
                             levels = sum_df$Description[order(sum_df$n)])

# Use the order set above to set order of factor levels in mapped data
map$Description <- factor(map$Description, 
                          levels = sum_df$Description[order(sum_df$n)])

# look at distribution of annotations ----------------------

bin_plots<- list()
for(i in 1:length(as.character(unique(map$bin)))){
  query <- as.character(unique(map$bin))[i]
  plot <- ggplot(map %>% filter(bin == query), aes(x = Description, fill = origin)) + 
                geom_bar() + theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                xlab("KEGG Ortholog BRITE Hierarchy") + 
                ggtitle(query) +
                coord_flip() 
  bin_plots[[i]] <- plot
}

pdf("plots.pdf", onefile = TRUE)
bin_plots
dev.off()

for(i in 1:length(as.character(unique(map$bin)))){
  query <- as.character(unique(map$bin))[i]
  print(query)
  print(map %>%
         filter(bin == query) %>%
         filter(origin == "crumb") %>%
         tally()) 
}

# transfer taxonomy information 
tax <- merge(info, map, by.x = "name", by.y = "bin")
# color by kingdom 
plot2b <- ggplot(tax, aes(x = Description, fill = kingdom)) + 
                  geom_bar() + theme_minimal() +
                  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                  xlab("KEGG Ortholog BRITE Hierarchy") + 
                  ggtitle("Total # Occurrences of any Ortholog") +
                  coord_flip() 

plot2b

# total number of KOs in plass assembly vs. KOs in bin, per bin ----------------
num_tot_bin <- function(kegg_df, Bin){
  hu <- kegg_df %>% filter(bin == Bin) # filter to bin of interest
  hu_bin <- hu %>% filter(origin == "bin") # filter to hu bins
  return(length(hu_bin$geneID)) # take all
}

num_tot_crumb <- function(kegg_df, Bin){
  hu <- kegg_df %>% filter(bin == Bin)
  hu_crumb <- hu %>% filter(origin == "crumb") # subset out crumb
  return(length(hu_crumb$geneID)) # take all
}

num_tot_bin1 <- sapply(unique(kegg_df$bin), function(x){num_tot_bin(kegg_df = kegg_df, Bin = x)})
num_tot_crumb1 <- sapply(unique(kegg_df$bin), function(x){num_tot_crumb(kegg_df = kegg_df, Bin = x)})
hist(num_tot_crumb1/num_tot_bin1)

# number of unique KOs in plass assembly vs. unique KOs in bin, per bin --------

# bin
num_unique_bin <- function(kegg_df, Bin){
  hu <- kegg_df %>% filter(bin == Bin) # filter to bin of interest
  hu_bin <- hu %>% filter(origin == "bin") # filter to hu bins
  len <- length(unique(hu_bin$geneID)) # take unique (i.e. remove dups), get num of unique
  return(len) # return num
}

# crumb
num_unique_crumb <- function(kegg_df, Bin){
  hu <- kegg_df %>% filter(bin == Bin)
  hu_crumb <- hu %>% filter(origin == "crumb") # subset out crumb
  hu_bin <- hu %>% filter(origin == "bin") # subset out bin
  df<-hu_crumb[!(hu_crumb$geneID %in% hu_bin$geneID), ] # get KOs in crumb that do not occur in bin
  return(length(unique(df$geneID))) # take "unique" of these KOs (i.e. remove dups), return num
}

num_bin <- sapply(unique(kegg_df$bin), function(x){num_unique_bin(kegg_df = kegg_df, Bin = x)})
num_crumb <- sapply(unique(kegg_df$bin), function(x){num_unique_crumb(kegg_df = kegg_df, Bin = x)})
hist(num_crumb/num_bin)
num_crumb
num_bin
# grab unique KOs per bin --------------------------------

# bin
unique_bin <- function(kegg_df, Bin){
  hu <- kegg_df %>% filter(bin == Bin) # filter to bin of interest
  hu_bin <- hu %>% filter(origin == "bin") # filter to hu bins
  uni <- unique(hu_bin$geneID) # take unique (i.e. remove dups), get num of unique
  df <- data.frame("geneID" = uni, "bin" = rep(Bin, length(uni)), "origin" = rep("bin", length(uni)))
  return(df)
}

# crumb
unique_crumb <- function(kegg_df, Bin){
  hu <- kegg_df %>% filter(bin == Bin)
  hu_crumb <- hu %>% filter(origin == "crumb") # subset out crumb
  hu_bin <- hu %>% filter(origin == "bin") # subset out bin
  uni <- !(hu_crumb$geneID %in% hu_bin$geneID) # get KOs in crumb that do not occur in bin
  uni <- unique(df$geneID) # take "unique" of these KOs (i.e. remove dups), return num
  df <- data.frame("geneID" = uni, "bin" = rep(Bin, length(uni)), "origin" = rep("plass", length(uni)))
  return(df)
  }

