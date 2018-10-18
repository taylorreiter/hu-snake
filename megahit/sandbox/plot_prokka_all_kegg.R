# plot functional categories from kegg

# read in data
all_kegg_full <- list.files("~/Desktop/hu-croissant-unitigs-prokka-all-kegg/", ".txt$", full.names = T)
all_kegg <- list.files("~/Desktop/hu-croissant-unitigs-prokka-all-kegg/", ".txt$")
all_kegg <- gsub(".ko-ann.txt", "", all_kegg)
kegg_df <- data.frame(locus_tag = NA, geneID = NA, name = NA)
for(i in 1:(length(all_kegg_full))){
  kegg <- read.csv(all_kegg_full[i], sep = "\t", col.names = c("locus_tag", "geneID"), 
                   stringsAsFactors = FALSE, na.strings = "") 
  kegg$name <- rep(all_kegg[i], nrow(kegg))
  kegg_df <- rbind(kegg_df, kegg)
}
kegg_df


# Get mappings of KO id to BRITE hierarchy description
library(clusterProfiler)
keggo <- kegg_df[ , 2]
keggo <- keggo[!is.na(keggo)]
# set p value to 1 so everything will be returned;
# i.e. this is not actually an enrichment analysis
kegg_enrich <- enrichKEGG(gene = keggo, organism = "ko", pvalueCutoff = 1) 
# grab data.frame of results
df <- kegg_enrich@result




# Plot the total occurence of each KO id. e.g. each KO id is only represented once. 
# Total number of samples containing at least one member
library(ggplot2)
library(ggthemes)
occur <- data.frame(id = rep(df$Description, df$Count))
ggplot(occur, aes(x = id)) + geom_bar() + theme_few() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("KEGG Ortholog BRITE Hierarchy") + 
  ggtitle("Total # of Hu Croissants Containing at least One Ortholog") +
  # scale_y_continuous(breaks = c(5, 10), labels = c("5", "10"))
  coord_flip()

sum_df <- occur %>%
  group_by(id) %>%
  tally() 

sum_df$id<- factor(sum_df$id, 
                             levels = sum_df$id[order(sum_df$n)])

ggplot(data = sum_df, aes(x = id, y = n)) +
  geom_bar(stat = "identity") + theme_few() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("KEGG Ortholog BRITE Hierarchy") +
  ylab("Count") +
  ggtitle("Functional Categories Contained in Hu Croissant Unitigs") +
  coord_flip()


# Plot total # of occurrences of any member in a sample and for each category
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
  colnames(kegg_df) <- c("locus_tag", "geneID")
  kegg_df <- right_join(kegg_df, mapping)
  kegg_df <- kegg_df[!is.na(kegg_df$locus_tag), ]
  kegg_df <- kegg_df[!is.na(kegg_df$geneID), ]
  head(kegg_df)
  kegg_df$sample <- rep(df_label, nrow(kegg_df))
  return(kegg_df)
}

map08 <- transfer_mappings(kegg_df = kegg08, mappings = mappings, df_label = "kegg08")
map11 <- transfer_mappings(kegg_df = kegg11, mappings = mappings, df_label = "kegg11")
map19 <- transfer_mappings(kegg_df = kegg19, mappings = mappings, df_label = "kegg19")

lab_df <- rbind(map08, map11, map19)
sum_df <- lab_df %>%
  group_by(Description) %>%
  tally() 

sum_df$Description <- factor(sum_df$Description, 
                             levels = sum_df$Description[order(-sum_df$n)])

# ggplot(data = sum_df, aes(x = Description, y = n, fill = sample)) + 
#   geom_bar(stat = "identity") + theme_few() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   xlab("KEGG Ortholog BRITE Hierarchy") +
#   ylab("Count") + 
#   ggtitle("Functional Categories Contained in Hu Croissants") +
#   scale_y_continuous(breaks = c(5, 10), labels = c("5", "10")) 

# Use the order set above to order old DF. 
lab_df$Description <- factor(lab_df$Description, 
                             levels = sum_df$Description[order(-sum_df$n)])
library(ggplot2)
library(ggthemes)
ggplot(lab_df, aes(x = Description, fill = sample)) + geom_bar() + theme_few() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("KEGG Ortholog BRITE Hierarchy") + 
  ggtitle("Total # Ocurrences of any Ortholog in a Functional Category in 3 Hu Croissants") +
  scale_y_continuous(breaks = c(5, 10), labels = c("5", "10")) 


