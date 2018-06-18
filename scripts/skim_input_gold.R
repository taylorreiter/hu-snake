library(dplyr)
library(tidyr)
library(skimr)

# list files
all_fai <- list.files("inputs/hu-crumbs_gold", ".fai$")
all_fai_full <- list.files("inputs/hu-crumbs_gold", ".fai$", full.names = T)

# read in all files
fai_df <- data.frame(V2 = NA, name = NA)
for(i in 1:(length(all_fai))){
  fai <- read.table(all_fai_full[i], head = F, stringsAsFactors = F)
  fai$name <- rep(all_fai[i], nrow(fai))
  fai <- fai[ , c("V2", "name")]
  fai_df <- rbind(fai_df, fai)
}

# change names of files
fai_df[] <- lapply(fai_df, function(x) gsub(".fa.cdbg_ids.contigs.fa.gz.crumbs.fa.fai", ";unitig", x))
fai_df[] <- lapply(fai_df, function(x) gsub(".fa.cdbg_ids.contigs.fa.gz.crumbs.fa.sub.fa.fai", ";sub", x))
fai_df[] <- lapply(fai_df, function(x) gsub(".fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa.fai", ";assembly", x))

# Remove NA
fai_df <- na.omit(fai_df)

# separate genome number from type using semi colon as delimiter
fai_df <- fai_df %>%
          separate(name, c("genome", "type"), ";")

colnames(fai_df) <- c("length", "genome", "type")

# convert to factor
fai_df$genome <- as.factor(fai_df$genome)
fai_df$type <- as.factor(fai_df$type)
fai_df$length <- as.numeric(fai_df$length)

# Change default skim stats so sum (aka number of nucleotides) will be included.
skim_with(numeric = list(sum = sum), append = TRUE)
print("skim!")
skim_summary <- fai_df %>% 
  group_by(type, genome) %>%
  skim_to_wide() %>%
  select(genome, n, p0, p50, p100, sum, hist)

write.table(skim_summary, file = "outputs/hu-crumbs_gold/summary_of_inputs.tsv", sep = "\t", quote = F, row.names = F)
