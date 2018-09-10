library(Biostrings)
library(tidyverse)
setwd('~/github/hu-snake')

sb1 <- readAAStringSet("outputs/hu-crumbs-bin/assembly/prokka/sb1.faa")
sb1 <- readDNAStringSet("outputs/hu-crumbs-bin/assembly/prokka/sb1.ffn")
sb1 <- readDNAStringSet("outputs/hu-crumbs-bin/assembly/megahit/sb1.contigs.fa") # 222 duplicates
sb1 <- readDNAStringSet("sandbox/ndj_uni/sb1.fa")
length(sb1)

sb1 <- sort(sb1)
ndj <- sb1[duplicated(sb1)]
ndj <- sort(ndj)

dups <- list()
for(i in 1:length(ndj)){
  dups[[i]] <- names(sb1[which(vcountPattern(ndj[[i]], sb1) == 1)])
}

 unlist(dups)
tmp <- Reduce(rbind, dups)
tmp <- as.data.frame(tmp)
tmp$bin1 <- str_extract(pattern = "[^\\_]*", tmp$V1)
tmp$bin1 <- gsub("ass", "", tmp$bin1)
tmp$bin2 <- str_extract(pattern = "[^\\_]*", tmp$V2)
tmp$bin2 <- gsub("ass", "", tmp$bin2)
tmp <- tmp %>%
  group_by(bin1, bin2) %>%
  tally()

map<- data.frame(bin = info$name, tax = info$taxonomy)

tmp <- left_join(tmp, map, by = c("bin1" = "bin"))
tmp <- left_join(tmp, map, by = c("bin2" = "bin"))
tmp
colnames(tmp) <- c("bin1", "bin2", "count", "tax1", "tax2")



# on aas ------------------------------------------------------------------

library(Biostrings)
library(tidyverse)
setwd('~/github/hu-snake')

sb12 <- readAAStringSet("outputs/hu-crumbs-bin/assembly/prokka/sb1.faa")
sb1 <- readDNAStringSet("outputs/hu-crumbs-bin/assembly/prokka/sb1.ffn")
sb1 <- readDNAStringSet("outputs/hu-crumbs-bin/assembly/megahit/sb1.contigs.fa") # 222 duplicates
length(sb12)

sb12 <- sort(sb12)
ndj2 <- sb12[duplicated(sb12)]
ndj2 <- sort(ndj2)

dups2 <- list()
for(i in 1:length(ndj2)){
  dups2[[i]] <- names(sb12[which(vcountPattern(ndj2[[i]], sb12) == 1)])
}

unlist(dups2)
tmp2 <- Reduce(rbind, dups2)
tmp2 <- as.data.frame(tmp2)
tmp2$bin1 <- str_extract(pattern = "[^\\_]*", tmp2$V1)
tmp2$bin1 <- gsub("ass", "", tmp2$bin1)
tmp2$bin2 <- str_extract(pattern = "[^\\_]*", tmp2$V2)
tmp2$bin2 <- gsub("ass", "", tmp2$bin2)
tmp2 <- tmp2 %>%
  group_by(bin1, bin2) %>%
  tally()


tmp2 <- left_join(tmp2, map, by = c("bin1" = "bin"))
tmp2 <- left_join(tmp2, map, by = c("bin2" = "bin"))
tmp2


