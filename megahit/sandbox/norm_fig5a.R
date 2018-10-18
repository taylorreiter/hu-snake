# the directory in which this R script lives contains sym-linked files of the crumb assembly GhostKOALA output and the bin GhostKOALA output. 
# the purpose of this script is to visualize the relative impact of the annotations from the crumbs against that of the bins.

library(clusterProfiler)
library(Biostrings)
library(dplyr)
library(ggplot2)

# GET CONTIG NAMES OF DUPLICATED NUCLEOTIDE SEQUENCES ---------------------------

sb1 <- readDNAStringSet("outputs/hu-crumbs-bin/assembly/prokka/sb1.ffn")
ndj <- sb1[duplicated(sb1)]

dups2 <- list()
for(i in 1:length(ndj)){
  dups2[[i]] <- names(sb1[which(vcountPattern(ndj[[i]], sb1) == 1)])
}

keep <- vector()
for(i in 1:length(dups2)){
  keep[i] <- dups2[[i]][1]
}

remove <- !(unlist(dups2) %in% keep) 
remove <- unlist(dups2)[remove]
remove <- gsub(" .*", "", remove)

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

info <- read.csv("inputs/hu_info.csv") # read in metadata
info <- info %>% 
         filter(sample_origin == "SB1") # only keep SB1 samples


crumb_df <- import_kegg(file = "outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt", 
                        gsub_regex = "(ass_[0-9]*)", origin = "crumb")

crumb_df <- crumb_df %>% 
              filter(bin %in% info$name) # filter to SB1 samples

bin_df <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", 
                      gsub_regex = "(_[0-9]*)", origin = "bin")
kegg_df <- rbind(crumb_df, bin_df) # combine two dfs

kegg_df <- kegg_df %>%
            filter(bin %in% crumb_df$bin) %>% # filter out the extra hu bins
            filter(!(locus_tag %in% remove)) # remove idential nucleotide (non disjoint)

keggo <- kegg_df[ , 2] # grab kegg orthologs
keggo <- keggo[!is.na(keggo)] # remove NAs

kegg_enrich <- enrichKEGG(gene = keggo, organism = "ko", pvalueCutoff = 1) # set p value to 1 so everything will be returned

df <- kegg_enrich@result # grab data.frame of results

# PLOT  -------------------------------------------------------

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
map <- transfer_mappings(kegg_df = kegg_df, mappings = mapping, df_label = "all")

sum_df <- map %>%
            group_by(Description) %>%
            tally() # tally the description; use to set order of the plot

sum_df <- sum_df[order(sum_df$n, decreasing = T), ] # order by total occurence
sum_df <- sum_df[1:20, ] # retain only top 20
#sum_df$Description <- droplevels(sum_df$Description)
sum_df$Description <- factor(sum_df$Description, 
                             levels = sum_df$Description[order(sum_df$n)]) # set order from overall sums

map <- map %>% 
  filter(Description %in% sum_df$Description) # prune map to only top 20

# Use the order set above to set order of factor levels in mapped data
map$Description <- factor(map$Description, 
                          levels = sum_df$Description[order(sum_df$n)])

bin <- map %>% 
          filter(origin == "bin") %>% 
          group_by(Description) %>%
          tally()

crumb <- map %>% 
          filter(origin == "crumb") %>% 
          group_by(Description) %>%
          tally()

all <- map %>%
        group_by(Description) %>% 
        tally()

bin$dist <- bin$n/all$n
bin$dist2 <- bin$n/sum(bin$n)
bin$origin <- rep('bin', nrow(bin))
crumb$dist <- crumb$n/all$n
crumb$dist2 <- crumb$n/sum(crumb$n)
crumb$origin <- rep('crumb', nrow(crumb))

tmp <- rbind(bin, crumb)
plot <- ggplot(tmp, aes(x = Description, y = dist2, color = origin)) + 
                geom_point() + theme_minimal() +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                xlab("ortholog pathways") + 
                ylab("# in pathway / total observations, for either crumb or bin") +
                ggtitle("Total # Occurrences of any Ortholog") +
                scale_fill_hue(labels = c("binned genomes", "unbinned content"))
plot


# make plot with two y axes
plot2 <- ggplot(tmp, aes(x = Description, y = dist, color = origin)) +
                geom_point() + theme_minimal() +
                theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                xlab("ortholog pathways") + 
                ggtitle("Total # Occurrences of any Ortholog") +
                scale_fill_hue(labels = c("binned genomes", "unbinned content")) +
                scale_y_continuous(sec.axis = sec_axis(trans = ~./10))
plot2           
