# the directory in which this R script lives contains sym-linked files of the crumb assembly GhostKOALA output and the bin GhostKOALA output. 
# the purpose of this script is to visualize the relative impact of the annotations from the crumbs against that of the bins.

setwd("/Users/taylorreiter/github/hu-snake/sandbox/bin_plus_crumb_assembly/")

library(clusterProfiler)
library(dplyr)
library(ggplot2)

# READ & FORMAT DATA ------------------------------------------------------------

# read in metadata
info <- read.csv("~/github/hu-snake/sandbox/hu_info.csv")

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

crumb_df <- import_kegg(file = "outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt", 
                        gsub_regex = "(ass_[0-9]*)", origin = "crumb")
bin_df <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", 
                      gsub_regex = "(_[0-9]*)", origin = "bin")

# combine two dfs
kegg_df <- rbind(crumb_df, bin_df)

# Get mappings of KO id to BRITE hierarchy description
keggo <- kegg_df[ , 2]
keggo <- keggo[!is.na(keggo)]

# To run for bin:
# keggo <- bin_df[ , 2]
# keggo <- keggo[!is.na(keggo)]

# to run for crumbs:
# keggo <- crumb_df[ , 2]
# keggo <- keggo[!is.na(keggo)]

# set p value to 1 so everything will be returned;
# i.e. this is not actually an enrichment analysis
kegg_enrich <- enrichKEGG(gene = keggo, organism = "ko", pvalueCutoff = 1) 
dotplot(kegg_enrich)
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

# write plot
#pdf(args[3], width=8, height=9)
plot1
#dev.off()

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
plot2a <- ggplot(map, aes(x = Description, fill = origin)) + geom_bar() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("KEGG Ortholog BRITE Hierarchy") + 
  ggtitle("Total # Ocurrences of any Ortholog") +
  coord_flip() +
  scale_fill_hue(labels = c("query", "crumb assembly"))

# write plot
#pdf(args[4], width = 10, height=9)
plot2a
#dev.off()


# look at distribution of crumb assembly annotations ----------------------

bin_plots<- list()
for(i in 1:length(as.character(unique(map$bin)))){
  query <- as.character(unique(map$bin))[i]
  plot <- ggplot(map %>% filter(bin == query), aes(x = Description, fill = origin)) + geom_bar() + theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    xlab("KEGG Ortholog BRITE Hierarchy") + 
    ggtitle(query) +
    coord_flip() +
    scale_fill_hue(labels = c("query", "crumb assembly")) 
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
#tax <- merge(info, map, by.x = "name", by.y = "bin")
# color by kingdom 
plot2b <- ggplot(tax, aes(x = Description, fill = kingdom)) + geom_bar() + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("KEGG Ortholog BRITE Hierarchy") + 
  ggtitle("Total # Occurrences of any Ortholog") +
  coord_flip() 

# write plot
#pdf(args[5], width=8, height=9)
plot2b
#dev.off()


# Explore! ----------------------------------------------------------------

# which archaea has peptidoglycan?
tax %>% filter(Description == 'Peptidoglycan biosynthesis')
# hu-genome30ass_00647
kegg_df %>% filter(locus_tag == 'hu-genome30ass_00647')  
# >hu-genome30ass_00647
# ATGCCCGAGGAGATCCACACGGACGATCCCCCCGATAGCAGGCGCAGGGGAAAGCACCGG
# GACAGGAGGTCAATCGGTCCGGTCCCGAACCTTGAAGAGCACGTGAAATCCGACTGGTGG
# AGAGGCATCTTCAACCGCCTCTACTTAAAGACCGACGTCGACGTCGTCGACGACCCGCGG
# ATCACCGAAAGAGAGATCGACCGGATCGCCCGGATCCTGCACCTCCAGCCCGATGAAAAG
# ATCCTTGACCTCTGCTGCGGCCAGGGAAGGCACACCCTCGAACTCGCGCGCAGGGGCTAC
# AACGCCGAAGGCCTCGACCAGTCGCACTACCTGATCCAGCGGGCCCGGTCGACCGCGAAG
# AAGGAGAGGCTCCCGGTCAGGTTCCGGGAGGGCGACGCACGCCGGCTCCCGTACCGCACC
# GACACCTGCGACGTCGTCCTCGTCCTCGGGAACAGTTTCGGCTACTTCGATTCCGTCGAG
# GAAGACCTGCGGATCCTCACCGAGCTCCGGCGTATCCTCAAAGCCGCGGAGGAGCGGTTC
# GGCATGGAGCCGGATACGCGGGAGACCGGCAAAGAACCTCTCCCGCTCGACCTTACCCGA
# CAGACGGCCTGA

# >hu-genome30ass_00647
# MPEEIHTDDPPDSRRRGKHRDRRSIGPVPNLEEHVKSDWWRGIFNRLYLKTDVDVVDDPR
# ITEREIDRIARILHLQPDEKILDLCCGQGRHTLELARRGYNAEGLDQSHYLIQRARSTAK
# KERLPVRFREGDARRLPYRTDTCDVVLVLGNSFGYFDSVEEDLRILTELRRILKAAEERF
# GMEPDTRETGKEPLPLDLTRQTA

# Match: methyltransferase domain-containing protein [Methanoculleus sp. MH98A] 

# We accept that there is some error in KEGG...?


# nif ---------------------------------------------------------------------

nifA <- 'K02584' #Nif-specific regulatory protein
nifB <- 'K02585' #nitrogen fixation protein NifB
nifD <- 'K02586' #nitrogenase molybdenum-iron protein alpha chain [EC:1.18.6.1]
nifE <- 'K02587' #nitrogenase molybdenum-cofactor synthesis protein NifE
nifH <- 'K02588' #nitrogenase iron protein NifH [EC:1.18.6.1]
nifHD1 <- 'K02589' #nifI1 nitrogen regulatory protein PII 1
nifHD2 <- 'K02590' #nifI2 nitrogen regulatory protein PII 2
nifK <- 'K02591' #nitrogenase molybdenum-iron protein beta chain [EC:1.18.6.1]
nifN <- 'K02592' #nitrogenase molybdenum-iron protein NifN
nifT <- 'K02593' #nitrogen fixation protein NifT
nifV <- 'K02594' #homocitrate synthase NifV [EC:2.3.3.14]
nifW <- 'K02595' #nitrogenase-stabilizing/protective protein
nifX <- 'K02596' #nitrogen fixation protein NifX
nifZ <- "K02597" # nitrogen fixation protein NifZ
glnB <- 'K04751' # nitrogen regulatory protein P-II 1
nifs <- c(nifA, nifB, nifD, nifE, nifH, nifHD1, nifHD2, nifK, nifN, nifT, nifV, nifW, nifX, nifZ, glnB)

kegg_df[kegg_df$geneID %in% nifs, ]


# marker genes ------------------------------------------------------------

gyrA <- 'K02469'
gyrB <- 'K02470'
recA <- 'K03553'

kegg_df[kegg_df$geneID %in% c(gyrA, gyrB, recA), ]
