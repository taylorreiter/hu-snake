# the directory in which this R script lives contains sym-linked files of the crumb assembly GhostKOALA output and the bin GhostKOALA output. 
# the purpose of this script is to visualize the relative impact of the annotations from the crumbs against that of the bins.

library(dplyr)
library(ggplot2)
library(tidyr)

# READ & FORMAT DATA ------------------------------------------------------------
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

info <- read.csv(snakemake@input[["info"]]) # read in metadata
#info <- read.csv("inputs/hu_info.csv")
info <- info %>% 
          filter(sample_origin == "SB1") # only keep SB1 samples

crumb_df <- import_kegg(file = snakemake@input[["crumb_kegg"]], 
                        gsub_regex = "(_SRR1976948.[0-9_]{2,12})", origin = "crumb")

# crumb_df <- import_kegg(file = "outputs/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_SRR1976948.[0-9_]{2,12})", origin = "crumb")
crumb_df <- crumb_df %>% 
              filter(bin %in% info$name) # filter to SB1 samples

bin_df <- import_kegg(snakemake@input[["bin_kegg"]], 
                      gsub_regex = "(_[0-9]*)", origin = "bin")
# bin_df <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_[0-9]*)", origin = "bin")
kegg_df <- rbind(crumb_df, bin_df) # combine two dfs

kegg_df <- kegg_df %>%
            filter(bin %in% crumb_df$bin) # filter out the extra hu bins

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
  uni <- hu_crumb[!(hu_crumb$geneID %in% hu_bin$geneID), ] # get KOs in crumb that do not occur in bin
  uni <- unique(uni$geneID) # take "unique" of these KOs (i.e. remove dups), return num
  df <- data.frame("geneID" = uni, "bin" = rep(Bin, length(uni)), "origin" = rep("crumb", length(uni)))
  return(df)
}

uni_crumb <- lapply(unique(info$name), function(x) unique_crumb(kegg_df = kegg_df, Bin = x))
uni_bin <- lapply(unique(info$name), function(x) unique_bin(kegg_df = kegg_df, Bin = x))
uni_crumb <- do.call(rbind, uni_crumb)
uni_bin <- do.call(rbind, uni_bin)

uni <- rbind(uni_crumb, uni_bin)
# DO KEGG -----------------------------------------------

keggparse <- read.delim("explore/ko00001_parse.txt", header= F, sep = "\t")
keggparse <- separate(keggparse, V4, into = c("geneID", "description"), sep = " ", remove= T, extra = "merge", fill = "left")
keggparse <- separate(keggparse, V3, into = c("path", "num"), sep = "\\[", remove = T, extra = "merge", file = "left")
keggparse <- separate(keggparse, path, into = c("num2", "path"), sep = " ", remove = T, extra = "merge", file = "left")
keggparse <- separate(keggparse, V2, into = c("num3", "category"), sep = " ", remove = T, extra = "merge", file = "left")

# deal with duplicates -- there are only 22,236 kegg orthologs in keggparse, but there are 49,935 rows, meaning there are duplicates. 

keggparse <- keggparse[match(unique(keggparse$geneID), keggparse$geneID),]

# keggo[!(keggo %in% keggparse$geneID)]
# levels(keggparse$V2) # 49
# levels(keggparse$V3) # 511

# PLOT 5a -------------------------------------------------------

# Transfer mappings to original kegg data
map <- left_join(uni, keggparse, by = "geneID")

sum_df <- map %>%
            group_by(path) %>%
            tally() # tally the description; use to set order of the plot

sum_df <- sum_df[order(sum_df$n, decreasing = T), ] # order by total occurence
sum_df <- sum_df[1:20, ] # retain only top 20
#sum_df$path <- droplevels(sum_df$path)
sum_df$path <- factor(sum_df$path, 
                             levels = sum_df$path[order(sum_df$n)]) # set order from overall sums

map <- filter(map, path %in% sum_df$path) # prune map to only top 20

# Use the order set above to set order of factor levels in mapped data
map$path <- factor(map$path, levels = sum_df$path[order(sum_df$n)])

plot5a <- ggplot(map, aes(x = path, fill = origin)) + geom_bar() + theme_minimal() +
                theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
                xlab("ortholog pathways") + 
                ggtitle("Unique Occurrences of any Ortholog") +
                coord_flip() +
                scale_fill_hue(labels = c("unbinned content", "binned genomes"))

# write plot
pdf(snakemake@output[["pdf"]], width = 6, height = 3)
#pdf('outputs/figures/fig5a.pdf', width = 6, height = 3)
plot5a
dev.off()

png(snakemake@output[["png"]], width = 6, height = 3, units = 'in', res = 300)
#png('outputs/figures/fig5a.png', width = 6, height = 3, units = 'in', res = 300)
plot5a
dev.off()

