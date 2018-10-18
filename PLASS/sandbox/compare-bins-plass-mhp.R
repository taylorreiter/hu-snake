# compare crumb assemblies to PLASS assemblies. 


# import data -------------------------------------------------------------

info <- read.csv("inputs/hu_info.csv")
info <- filter(info, sample_origin == "SB1")

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
plass <- filter(plass, bin %in% info$name)

mhp <- import_kegg(file = "../../outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt", gsub_regex = "(ass_[0-9]*)", origin = "mhp")
mhp <- filter(mhp, bin %in% info$name)

bin <- import_kegg(file = "outputs/hu-bins/GhostKOALA/user_ko_definition.txt", gsub_regex = "(_[0-9]*)", origin = "bin")
bin <- filter(bin, bin %in% info$name)




# sb1 level ---------------------------------------------------------------
# grab only unique kegg orthologs
plassg <- unique(plass$geneID)
mhpg <- unique(mhp$geneID)
bing <- unique(bin$geneID)

length(plassg[!(plassg %in% bing)]) 


length(plassg) - length(intersect(bing, plassg))
length(mhpg) - length(intersect(mhpg, plassg))
length(plassg) - length(intersect(plassg, mhpg))
length(bing) - length(intersect(bing, plassg))
length(bing) - length(intersect(bing, mhpg))

# in plass, not in bins: 136
# in plass, not in megahit/prokka: 1466
# in bins, not in plass: 7
# in bins, not in megahit/prokka: 1344


# per bin basis -----------------------------------------------------------

# see plot_per_bin_kegg.R
