library(KEGGREST)
library(tidyr)

# Parse An et al. genes

an <- read.csv("sandbox/an-et-al/an_et_al_genes.csv", stringsAsFactors = F)

# look up kegg orthlogs
ko <- lapply(an$EC, keggFind, database = 'ko')

# parse kegg id
ko <- unlist(ko)
kegg_ids <- names(ko)
kegg_ids <- gsub("ko:", "", kegg_ids)

# make df of id, name, and ec
kegg_df <- data.frame(ko = names(ko), gene = ko)
kegg_df <- separate(data = kegg_df, col = gene, into = c("symbol", "gene"), sep = ";")
kegg_df <- separate(data = kegg_df, col = gene, into = c("gene", "EC"), sep = "\\[")
kegg_df[] <- lapply(kegg_df, function(x) gsub("\\]", "", x))
kegg_df[] <- lapply(kegg_df, function(x) gsub("ko:", "", x))
kegg_df <- separate(data = kegg_df, col = EC, into = c("EC", "EC2", "EC3"), sep = " ", extra = "merge", fill = "right")

# Create single info dataframe
info <- merge(kegg_df, an)
# Check that only PFAM ids fail
an[!(an$EC %in% info$EC), ]

head(info)
info <- info[ , c("EC", "ko", "symbol", "gene", "EC2", "Enzyme", "type", "class", "oxygen")]
write.table(info, "sandbox/an-et-al/an_et_al_genes_parsed.tsv", quote = FALSE, row.names = FALSE, sep = "\t")
