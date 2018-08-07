# designate metabolic vs informational genes
# figure out how to randomly and proportionately assign genes to crumbs & queries
# "Calculate bias toward metabolic genes in crumbs" 
setwd("/Users/taylorreiter/github/hu-snake")
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
brite <- read.delim("sandbox/ko00001_parse.txt", header= F, sep = "\t")
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
crumb_full <- list.files("outputs/hu-crumbs_bin/assembly/GhostKOALA/", "full.txt$", full.names = T)
crumb_kegg <- list.files("outputs/hu-crumbs_bin/assembly/GhostKOALA/", "full.txt$")

crumb_kegg <- gsub(".ko-ann-full.txt", "", crumb_kegg)
crumb_df <- data.frame(locus_tag = NA, geneID = NA, name = NA, score = NA, second = NA, second_score = NA, bin = NA)
for(i in 1:(length(crumb_full))){
  kegg <- read.csv(crumb_full[i], sep = "\t", stringsAsFactors = FALSE, na.strings = "",
                   col.names = c("locus_tag", "geneID", "name", "score", "second", "second_score", "bin")) 
  crumb_df <- rbind(crumb_df, kegg)
}
# add row specifying identity
crumb_df$origin <- rep("crumb", nrow(crumb_df))

bin_full <- list.files("sandbox/bin_plus_crumb_assembly/bin", "full.txt$", full.names = T)
bin_full <- bin_full[2:31]
bin_kegg <- list.files("sandbox/bin_plus_crumb_assembly/bin/", "full.txt$")
bin_kegg <- bin_kegg[2:31]

bin_kegg <- gsub(".ko-ann-full.txt", "", bin_kegg)
bin_df <- data.frame(locus_tag = NA, geneID = NA, name = NA, score = NA, second = NA, second_score = NA, bin = NA)
for(i in 1:(length(bin_full))){
  kegg <- read.csv(bin_full[i], sep = "\t", stringsAsFactors = FALSE, na.strings = "",
                   col.names = c("locus_tag", "geneID", "name", "score", "second", "second_score", "bin")) 
  bin_df <- rbind(bin_df, kegg)
}

bin_df$origin <- rep("bin", nrow(bin_df))

kegg_df <- rbind(crumb_df, bin_df)

# retain only kegg ids for which the score is > 100
kegg_df <- kegg_df[kegg_df$score > 100, ]

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

# strawman part -----------------------------------------------------------

# 0 = informational
# 1 = metabolic
grep("Transcription", as.character(levels(map$Description)))
levels(droplevels(map$Description))
sample(c(0,1), size = nrow(map), replace = TRUE, prob = )
?sample
data = matrix(data = 1:6, ncol = 2, nrow = 6)
newData = cbind(data, sample(c(0,1), size = nrow(data), replace = TRUE))
