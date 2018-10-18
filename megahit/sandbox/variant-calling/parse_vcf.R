# Determine if the number of variants per bin correlates with crumbs size
# 1. Map gi accession to bin name. 
# 2. Read vcf file and merge header information
# 3. Count number of variants per bin

# 1. Map gi accession to bin name. 
# Get bin names for each bin
# cd ~/github/hu-snake/inputs/hu-genomes
# for infile in *fa; do grep ">" $infile > ${infile}.headers.txt; done

files <- list.files("../../inputs/hu-genomes", "headers.txt$", full.names = TRUE)
headers_all <- data.frame(gi = NA, V2 = NA, V3 = NA, V4 = NA, V5 = NA,  V6 = NA, V7 = NA, V8 = NA, V9 = NA, V10 = NA, V11 = NA, bin = NA)
for(infile in files){
  bin <- gsub(".fa.headers.txt", "", infile)
  bin <- gsub("\\.\\.\\/\\.\\.\\/inputs\\/hu-genomes\\/", "", bin)
  headers <- read.delim(infile, header = F, sep = " ", col.names = c("gi", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11"))
  headers$bin <- rep(bin, nrow(headers))
  headers_all <- rbind(headers_all, headers)
}

headers_all <- headers_all[ , c("gi", "bin")]
headers_all$gi <- sub(">", "", headers_all$gi)

# 2. Read in vcf file and merge with header info
vcf <- read.delim(file = "sb1_final_variants.vcf", header = F, comment.char = "#", stringsAsFactors = F)
vcf <- merge(headers_all, vcf, by.x = "gi", by.y = "V1")
head(vcf)
table(vcf$bin)

library(dplyr)
vcf_counts <- vcf %>% 
  group_by(bin) %>%
  tally()

# 3. Correlate with crumb metrics. 
#   a. Nucleotides in unitigs
#   b. crumb file size
#   c. 31-mers in crumbs

donut_files <- list.files("../../inputs/hu-s1_k31_r1_search_oh0_may20/", ".fa.cdbg_ids.contigs.fa.gz.donut.fa$", full.names = TRUE)

donut_bytes <- sapply(donut_files, file.size)
donut_bytes <- as.data.frame(donut_bytes)
donut_bytes$file <- row.names(donut_bytes)
donut_bytes$file <- gsub("\\.\\.\\/\\.\\.\\/inputs\\/hu-s1_k31_r1_search_oh0_may20\\/\\/", "", donut_bytes$file)
donut_bytes$file <- gsub(".fa.cdbg_ids.contigs.fa.gz.donut.fa", "", donut_bytes$file)


cor <- merge(vcf_counts, donut_bytes, by.x = "bin", by.y = "file")
plot(cor$n, cor$donut_bytes)

donut_unitig_counts <- vector()
for(infile in donut_files){
  count <- length(grep(">", readLines(infile)))
  donut_unitig_counts <- c(donut_unitig_counts, count)
  }

cor$donut_unitig_counts <- donut_unitig_counts
plot(cor$n, cor$donut_unitig_counts, xlab = "# Variants in Hu Bin", ylab = "# Unitigs in Donuts")

library(plotly)
plot_ly(cor, x = ~cor$n, y = ~cor$donut_unitig_counts, text = cor$bin) %>% 
  #add_markers(cor, color = ~cor$bin) %>% 
  layout(xaxis = list(title = "# Variants in Hu Bin"), yaxis = list(title = "# Unitigs in Donuts"))


# control for sequencing depth --------------------------------------------

inf <- read.csv("../hu_info.csv", stringsAsFactors = F)
head(inf)
norm <- merge(inf, vcf_counts, by.x = "name", by.y = "bin")
norm$counts_norm <- norm$n / norm$coverage 
norm <- merge(cor, norm, by.x = "bin", by.y = "name")
sb1 <- filter(norm, sample_origin == "SB1")
sb1
plot_ly(norm, x = ~sb1$n.x, y = ~sb1$coverage, text = sb1$bin) %>% 
  #add_markers(cor, color = ~cor$bin) %>% 
  layout(xaxis = list(title = "# Variants in Hu Bin "), yaxis = list(title = "coverage"))


# read in vcf correctly ---------------------------------------------------

library(vcfR)
vcfr <- read.vcfR("sb1_final_variants_trim_4-2.vcf")
chrom <- create.chromR(vcfr, name = "4:2")
chromoqc(chrom, dp.alpha=20)
