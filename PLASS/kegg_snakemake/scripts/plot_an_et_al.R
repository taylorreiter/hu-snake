#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# args[1]: INPUT. info -- "sandbox/an-et-al/an_et_al_genes_parsed.tsv"
# args[2]: INPUT. kegg directory -- "outputs/hu-croissant/*/GhostKOALA"; files end with *.ko-ann-full.txt
# args[3]: OUTPUT. tsv an_table --
# args[4]: OUTPUT. title of plot, as string 
# args[5]: OUTPUT. plot file name (e.g. "Subtracts")

library(dplyr)
library(ggplot2)
library(tidyr)

# check kegg ids for presence of enzymes indicated as important by An et al. 2013

# READ IN DATA ------------------------------------------------------------
# read in csv with EC numbers & KOs
print("read in info")
info <- read.delim("inputs/an_et_al_genes_parsed.tsv", stringsAsFactors = F, sep = "\t")

# function ----------------------------------------------------------------
print("defining functions")
find_an <- function(dir, kegg_ids, cutoff){
  # dir: directory where GhostKOALA output files are
  # kegg_ids
  # cutoff: Int. GhostKOALA cutoff score to use. Anything below this score will be removed.
  print("function: find_an()")
  all_kegg_anno_full <- list.files(dir, "full.txt$", full.names = T)
  all_kegg_anno <- list.files(dir, "full.txt$")
  all_kegg_anno <- gsub(".ko-ann-full.txt", "", all_kegg_anno)
  kegg_anno_df <- data.frame(V1 = NA, V2 = NA, V3 = NA, V4 = NA, name = NA)
  
  for(i in 1:(length(all_kegg_anno))){
    kegg_anno <- read.delim(all_kegg_anno_full[i], head = F, stringsAsFactors = F, quote = "")
    kegg_anno$name <- rep(all_kegg_anno[i], nrow(kegg_anno))
    kegg_anno <- kegg_anno[ , c("V1", "V2", "V3", "V4", "name")]
    kegg_anno_df <- rbind(kegg_anno_df, kegg_anno)
  }
  
  an <- kegg_anno_df[kegg_anno_df$V2 %in% kegg_ids, ]
  an <- an[an$V4 >= cutoff, ]
  return(an)
}

table_an <- function(an, info){
  # an: output of table_an function
  # info: information dataframe made at beginning of script. Must contain same column names used below
  print("function: table_an()")
  an <- merge(an, info, by.x = "V2", by.y = "ko")
  an <- an[ , c("V2", "name", "EC", "Enzyme", "symbol", "type", "class", "oxygen")]
  an <- an %>% 
    group_by(name, symbol, Enzyme, type, class, oxygen) %>%
    tally()
  an <- spread(an, key = name, value = n)
  return(an)
}

plot_an <- function(an, info, title){
  # an: output of table_an function
  # info: information dataframe made at beginning of script. Must contain same column names used below
  print("function: plot_an()")
  an <- merge(an, info, by.x = "V2", by.y = "ko")
  an <- an[ , c("V2", "name", "EC", "Enzyme", "symbol", "type", "class", "oxygen")]
  an <- an %>% 
    group_by(class, type, oxygen) %>%
    tally()
  plotty <- ggplot(an, aes(x = class, y = n, fill = type)) + 
    geom_bar(stat = "identity") + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylab("count") + ggtitle(title)
  return(plotty)
}

# apply -------------------------------------------------------------------
print("apply functions")
an_table <- table_an(an = find_an(dir = args[2], kegg_ids = info$ko, cutoff = 40), 
                 info = info) 

write.table(an_table, file = args[3], quote = FALSE, sep = "\t", row.names = F)

an <- find_an(dir = args[2], kegg_ids = info$ko, cutoff = 40)

an_plot <- plot_an(an = an, info = info, title = args[4])
# write plot
print("writing plot!")
pdf(file = args[5], width = 6, height = 5)
an_plot
dev.off()
