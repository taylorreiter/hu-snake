# check kegg ids for presence of enzymes indicated as important by An et al. 2013

library(dplyr)
library(ggplot2)

info <- read.csv("sandbox/an-et-al/an_et_al_genes_parsed.csv")
# function ----------------------------------------------------------------
find_an <- function(dir, kegg_ids, cutoff){
  # dir: directory where GhostKOALA output files are
  # kegg_ids
  # cutoff: Int. GhostKOALA cutoff score to use. Anything below this score will be removed.
  
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



subs <- table_an(an = find_an(dir = "~/Desktop/hu-croissant-subtract-prokka-all-kegg/", 
                              kegg_ids = kegg_ids, cutoff = 40), 
                 info = info)
uni <- table_an(an = find_an(dir = "~/Desktop/hu-croissant-unitigs-prokka-all-kegg/", 
                             kegg_ids = kegg_ids, cutoff = 40), 
                info = info)
assembly <- table_an(an = find_an(dir = "~/Desktop/hu-croissant-assemblies-prokka-kegg/", 
                                  kegg_ids = kegg_ids, cutoff = 40), 
                     info = info)
write.table(assembly, file = "sandbox/assembly_content_an2013.tsv", quote = FALSE, sep = "\t", row.names = F)
write.table(subs, file = "sandbox/subtracts_content_an2013.tsv", quote = FALSE, sep = "\t", row.names = F)
write.table(uni, file = "sandbox/unitig_content_an2013.tsv", quote = FALSE, sep = "\t", row.names = F)



subs <- find_an(dir = "~/Desktop/hu-croissant-subtract-prokka-all-kegg/", kegg_ids = kegg_ids, cutoff = 40)
uni <- find_an(dir = "~/Desktop/hu-croissant-unitigs-prokka-all-kegg/", kegg_ids = kegg_ids, cutoff = 40)
assembly <- find_an(dir = "~/Desktop/hu-croissant-assemblies-prokka-kegg/", kegg_ids = kegg_ids, cutoff = 40)


plot_an(an = uni, info = info, title = "Unitigs")
plot_an(an = sub, info = info, title = "Subtracts")
plot_an(an = assembly, info = info, title = "Assemblies")
