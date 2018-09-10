# designate metabolic vs informational genes
# figure out how to randomly and proportionately assign genes to crumbs & queries
# "Calculate bias toward metabolic genes in crumbs" 
setwd("/Users/taylorreiter/github/hu-snake/sandbox/PLASS/")
library(clusterProfiler)
library(tidyverse)

# get kegg hierarchy ------------------------------------------------------

# http://merenlab.org/2018/01/17/importing-ghostkoala-annotations/
# download.file('https://www.genome.jp/kegg-bin/download_htext?htext=ko00001.keg&format=htext&filedir=', destfile = "sandbox/ko00001.keg")
# in bash:
# kegfile="sandbox/ko00001.keg"
# 
# while read -r prefix content
# do
#   case "$prefix" in A) col1="$content";; \
#                     B) col2="$content" ;; \
#                     C) col3="$content";; \
#                     D) echo -e "$col1\t$col2\t$col3\t$content";;
#   esac 
# done < <(sed '/^[#!+]/d;s/<[^>]*>//g;s/^./& /' < "$kegfile") > sandbox/KO_Orthology_ko00001.txt
# There was a problem with special characters where the table was truncated. I opened it in excel and exported as a tsv and it solved these problems. 

# Read in ko table
brite <- read.delim("ko00001_parse.txt", header= F, sep = "\t")
brite <- separate(brite, V4, into = c("id", "description"), sep = " ", remov= T, extra = "merge", fill = "left")

gi <- brite %>% 
        filter(V1 == "09120 Genetic Information Processing")

bi <- brite %>%
        filter(V2 %in% c("09182 Protein families: genetic information processing", "09121 Transcription", 
                         "09122 Translation", "09123 Folding, sorting and degradation", "09124 Replication and repair", 
                         "09192 Unclassified: genetic information processing"))#, "09184 RNA family"))
info <- rbind(gi, bi)
# make sure each ko is only represented once 
info <- info[match(unique(info$id), info$id),]

# Read in all bin and crumb annotations
hu_info <- read.csv("inputs/hu_info.csv")
hu_info <- hu_info %>% 
            filter(sample_origin == "SB1")

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

plass <- import_kegg(file = "outputs/GhostKOALA/user_ko_definition.txt", 
                     gsub_regex = "(_SRR1976948.[0-9_]{2,12})", origin = "plass")
plass <- plass %>% 
          filter(bin %in% hu_info$name)


bin <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_[0-9]*)", origin = "bin")
bin <- hubin %>% 
          filter(bin %in% hu_info$name)

kegg_df <- rbind(plass, bin)

# merge info
kegg_df <- left_join(x = kegg_df, y = info, by = c("geneID" = "id"))

# remove locus_tags for which there is no kegg annotation
kegg_df <- kegg_df[!is.na(kegg_df$geneID), ]

# add a column to designate "info" and "other" for type of annotation. 
kegg_df$type <- ifelse(is.na(kegg_df$V1) == FALSE, "info", "other")

kegg_df %>% 
  group_by(origin, type) %>% 
  tally 

tmp <- kegg_df %>% 
  filter(origin == "crumb") %>%
  filter(type == "info")
33.7

# strawman part -----------------------------------------------------------

# 0 = informational
# 1 = metabolic
grep("Transcription", as.character(levels(map$Description)))
levels(droplevels(map$Description))
sample(c(0,1), size = nrow(map), replace = TRUE, prob = )
?sample
data = matrix(data = 1:6, ncol = 2, nrow = 6)
newData = cbind(data, sample(c(0,1), size = nrow(data), replace = TRUE))
