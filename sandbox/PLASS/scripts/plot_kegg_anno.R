#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# args[1]: INPUT. info -- "sandbox/hu_info.csv"
# args[2]: INPUT. kegg directory -- "outputs/hu-croissant/*/GhostKOALA"; files end with *.ko-ann-full.txt
# args[3]: OUTPUT. plot1
# args[4]: OUTPUT. plot2a
# args[5]: OUTPUT. plot2b  

library(clusterProfiler)
library(dplyr)
library(ggplot2)

# READ & FORMAT DATA ------------------------------------------------------------

# read in metadata
info <- read.csv(args[1])

# read in data
# the first sys arg will be a parameter which specifies in which folder the files are in. 
all_kegg_full <- list.files(args[2], "full.txt$", full.names = T)
all_kegg <- list.files(args[2], "full.txt$")

all_kegg <- gsub(".ko-ann-full.txt", "", all_kegg)
kegg_df <- data.frame(locus_tag = NA, geneID = NA, name = NA, score = NA, second = NA, second_score = NA, bin = NA)
for(i in 1:(length(all_kegg_full))){
  kegg <- read.csv(all_kegg_full[i], sep = "\t", stringsAsFactors = FALSE, na.strings = "",
                   col.names = c("locus_tag", "geneID", "name", "score", "second", "second_score", "bin")) 
  kegg$name <- rep(all_kegg[i], nrow(kegg))
  kegg_df <- rbind(kegg_df, kegg)
}

# retain only kegg ids for which the score is > 40
kegg_df <- kegg_df[kegg_df$score > 40, ]
print("Pruned ko ids at 40")

# Get mappings of KO id to BRITE hierarchy description
keggo <- kegg_df[ , 2]
keggo <- keggo[!is.na(keggo)]

# set p value to 1 so everything will be returned;
# i.e. this is not actually an enrichment analysis
kegg_enrich <- enrichKEGG(gene = keggo, organism = "ko", pvalueCutoff = 1) 
# grab data.frame of results
df <- kegg_enrich@result
print("obtained BRITE hierarchies")


# PLOT 1 ------------------------------------------------------------------
# Plot the total occurrence of each KO id. e.g. each KO id is only represented once. 

# Total number of samples containing at least one member
occur <- data.frame(id = rep(df$Description, df$Count))

# make a dataframe that counts each BRITE hierarchy
sum_df <- occur %>%
  group_by(id) %>%
  tally() 
print("plot 1 BRITE hierarchies tallied")

# set factor level to order plot by number of occurences
sum_df$id<- factor(sum_df$id, 
                   levels = sum_df$id[order(sum_df$n)])

# Plot
plot1 <- ggplot(data = sum_df, aes(x = id, y = n)) +
  geom_bar(stat = "identity") + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("KEGG Ortholog BRITE Hierarchy") +
  ylab("Count") +
  ggtitle("Total occurrence of each KEGG Ortholog") +
  coord_flip()

print("Writing plot1")

# write plot
pdf(args[3], width=8, height=9)
plot1
dev.off()

# PLOT 2 ------------------------------------------------------------------
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
  kegg_df <- right_join(kegg_df, mappings)
  kegg_df <- kegg_df[!is.na(kegg_df$locus_tag), ]
  kegg_df <- kegg_df[!is.na(kegg_df$geneID), ]
  kegg_df$sample <- rep(df_label, nrow(kegg_df))
  return(kegg_df)
}


# Transfer mappings
print("transfer mappings")
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

# color by genome bin (basically unintelligable)
plot2a <- ggplot(map, aes(x = Description, fill = name)) + geom_bar() + theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            xlab("KEGG Ortholog BRITE Hierarchy") + 
            ggtitle("Total # Ocurrences of any Ortholog") +
            coord_flip() 
plot2a
# write plot
pdf(args[4], width=10, height=9)
plot2a
dev.off()

# transfer taxonomy information 
tax <- merge(info, map, by = "name")
# color by kingdom 
plot2b <- ggplot(tax, aes(x = Description, fill = kingdom)) + geom_bar() + theme_minimal() +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            xlab("KEGG Ortholog BRITE Hierarchy") + 
            ggtitle("Total # Occurrences of any Ortholog") +
            coord_flip() 
# write plot
pdf(args[5], width=8, height=9)
plot2b
dev.off()

