library(clusterProfiler)
library(dplyr)
library(Biostrings)


# GET DUPLICATED NUC SEQS -------------------------------------------------

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

info <- read.csv(snakemake@input[["info"]]) # read in metadata
#info <- read.csv("inputs/hu_info.csv")
info <- info %>% 
          filter(sample_origin == "SB1") # only keep SB1 samples

crumb_df <- import_kegg(file = snakemake@input[["crumb_kegg"]], 
                        gsub_regex = "(ass_[0-9]*)", origin = "crumb")
#crumb_df <- import_kegg(file = "outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt", gsub_regex = "(ass_[0-9]*)", origin = "crumb")

crumb_df <- crumb_df %>% 
              filter(bin %in% info$name) %>% # filter to SB1 samples
              filter(!(locus_tag %in% remove)) # remove idential nucleotide (non disjoint)

keggo <- crumb_df[ , 2] # grab kegg orthologs
keggo <- keggo[!is.na(keggo)] # remove NAs
kegg_enrich <- enrichKEGG(gene = keggo, organism = "ko") # enrichment analysis
plot5b <- dotplot(kegg_enrich) # plot enrichment

pdf(snakemake@output[["pdf"]], width = 7.5, height = 3.5)
#pdf('outputs/figures/fig5b.pdf', width = 7.5, height = 3.5)
plot5b
dev.off()

png(snakemake@output[["png"]], width = 7.5, height = 3.5, units = 'in', res = 300)
#png('outputs/figures/fig5b.png', width = 7.5, height = 3.5, units = 'in', res = 300)
plot5b
dev.off()
