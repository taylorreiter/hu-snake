library(dplyr)
library(cluster)
library(dendextend)
library(factoextra)

comp <- read.csv(snakemake@input[['comp']])
rownames(comp) <- colnames(comp)
mat <- as.matrix(comp)
dend <- as.dendrogram(hclust(dist(mat), method = "ave"))
labels(dend) <- gsub("outputs.bam_subsets.", "", labels(dend))
labels(dend) <- gsub(".fa", "", labels(dend))

groups <- cutree(dend, h = 2)

# plot dendogram
#dend1 <- color_branches(dend, k = k)
#dend1 <- color_labels(dend1, k = k)
#plot(dend1, main = "default cuts by cutree")

for(i in unique(groups)){
  group <- names(groups[groups == i])
  write.table(group, file = paste0("outputs/comp/cut", i, ".txt"), 
              quote = F, row.names = F, col.names = F)
}
