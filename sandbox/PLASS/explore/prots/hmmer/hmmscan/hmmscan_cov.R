setwd("github/hu-snake/sandbox/PLASS/")
library(rhmmer)
library(GenomicRanges)
library(ggbio)
library(dplyr)

info <- read.csv("inputs/hu_info.csv")
info <- info %>% 
          filter(sample_origin == "SB1")

gyraplass <- read_domtblout('explore/prots/hmmer/plass-gyra-dom-hmmscan.out')
gyrabin <- read_domtblout('explore/prots/hmmer/bin-gyra-dom-hmmscan.out')
gyrbplass <- read_domtblout('explore/prots/hmmer/plass-gyrb-dom-hmmscan.out')
gyrbbin <- read_domtblout('explore/prots/hmmer/bin-gyrb-dom-hmmscan.out')

gyraplass <- gyraplass %>%
              mutate(bin = gsub("(_SRR1976948.[0-9_]{2,12})", "", query_name)) %>% 
              mutate(origin = rep("plass", nrow(gyraplass))) %>% 
              mutate(gene = rep("gyrA", nrow(gyraplass))) %>%
              filter(domain_ievalue < 1e-3) # suggested cutoff in hmmer manual

gyrabin <- gyrabin %>%
              mutate(bin = gsub("(_[0-9]*)", "", query_name)) %>% 
              mutate(origin = rep("bin", nrow(gyrabin))) %>% 
              mutate(gene = rep("gyrA", nrow(gyrabin))) %>%
              filter(domain_ievalue < 1e-3) %>%
              filter(bin %in% info$name)

markersA <- rbind(gyraplass, gyrabin) 

gyrbplass <- gyrbplass %>%
              mutate(bin = gsub("(_SRR1976948.[0-9_]{2,12})", "", query_name)) %>% 
              mutate(origin = rep("plass", nrow(gyrbplass))) %>% 
              mutate(gene = rep("gyrB", nrow(gyrbplass))) %>%
              filter(domain_ievalue < 1e-3)

gyrbbin <- gyrbbin %>%
            mutate(bin = gsub("(_[0-9]*)", "", query_name)) %>% 
            mutate(origin = rep("bin", nrow(gyrbbin))) %>% 
            mutate(gene = rep("gyrB", nrow(gyrbbin))) %>%
            filter(domain_ievalue < 1e-3) %>%
            filter(bin %in% info$name)

markersB <- rbind(gyrbplass, gyrbbin)

# plot
markersA <- makeGRangesFromDataFrame(markersA,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqinfo=NULL,
                                     seqnames.field=c("domain_name"),
                                     start.field="hmm_from",
                                     end.field="hmm_to",
                                     #strand.field="strand",
                                     starts.in.df.are.0based=FALSE)


ggplot(markers) + 
  stat_coverage(geom = "bar", facets = origin ~ seqnames) +
  theme_minimal() #+
  #scale_y_log10() +
  #ylab("log10 coverage")


plass_cov <- sapply(coverage(markers[markers$origin == "plass"]), as.numeric)
plass_cov <- data.frame(Chr=rep("gyrA", nrow(plass_cov)), Origin=rep("plass", nrow(plass_cov)), Position = 1:nrow(plass_cov), Coverage=as.numeric(unlist(plass_cov)))


bin_cov <- sapply(coverage(markers[markers$origin == "bin"]), as.numeric)
bin_cov <- data.frame(Chr=rep("gyrA", nrow(bin_cov)), Origin=rep("bin", nrow(bin_cov)), Position = 1:nrow(bin_cov), Coverage=as.numeric(unlist(bin_cov)))

covdf <- rbind(plass_cov, bin_cov)
ggplot(covdf, aes(Position, Coverage, fill=Origin)) + 
      geom_bar(stat="identity", position="identity") + facet_wrap(~Chr) +
      theme_minimal() +
      xlim(c(1, 583)) +
      scale_y_log10()


# calculate per bin gyrA -------------------------------------------------------

covdf_gyrA <- data.frame(pfam = character, Origin = character(), Position = numeric(), Coverage = numeric(), Bin = character())
for(hubin in unique(markersA$bin)){
    markers2 <- markers[markersA$bin == hubin]
    plass_cov <- sapply(coverage(markers2[markers2$origin == "plass"]), as.numeric)
    plass_cov <- data.frame(pfam=rep("gyrA", nrow(plass_cov)), 
                            Origin=rep("plass", nrow(plass_cov)), 
                            Position = 1:nrow(plass_cov), 
                            Coverage=as.numeric(unlist(plass_cov)), 
                            Bin = rep(hubin, nrow(plass_cov))) 
  
  bin_cov <- sapply(coverage(markers2[markers2$origin == "bin"]), as.numeric)
  if(typeof(nrow(bin_cov)) == "NULL"){
    covdf <- rbind(covdf, plass_cov)
    next
  } else {
  bin_cov <- data.frame(pfam=rep("gyrA", nrow(bin_cov)), 
                        Origin=rep("bin", nrow(bin_cov)), 
                        Position = 1:nrow(bin_cov), 
                        Coverage=as.numeric(unlist(bin_cov)), 
                        Bin = rep(hubin, nrow(bin_cov)))
  covdf_gyrA <- rbind(covdf_gyrA, plass_cov, bin_cov)
  }
}

covdf2_gyrA <- covdf_gyrA %>% filter(Coverage != 0)

ggplot(covdf_gyrA %>% filter(Coverage != 0), aes(Position, Coverage, fill=Origin)) + 
        geom_bar(stat="identity", position="identity") + facet_wrap(~ Bin) +
        theme_minimal() +
        xlim(c(1, 583)) + 
        scale_y_sqrt() + 
        ggtitle("Per Amino Acid Coverage of PFAM gyrA") +
        ylab("Square Root of Coverage Depth")


# calculate per bin gyrB --------------------------------------------------

markersB <- makeGRangesFromDataFrame(markersB,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqinfo=NULL,
                                     seqnames.field=c("domain_name"),
                                     start.field="hmm_from",
                                     end.field="hmm_to",
                                     #strand.field="strand",
                                     starts.in.df.are.0based=FALSE)

covdf_gyrB <- data.frame(pfam = character, 
                         Origin = character(), 
                         Position = numeric(), 
                         Coverage = numeric(), 
                         Bin = character())

for(hubin in unique(markersB$bin)){
  markers2 <- markersB[markersB$bin == hubin]
  plass_cov <- sapply(coverage(markers2[markers2$origin == "plass"]), as.numeric)
  if(typeof(nrow(plass_cov)) == "NULL"){
    next
  } else {
    plass_cov <- data.frame(pfam=rep("gyrB", nrow(plass_cov)), 
                            Origin=rep("plass", nrow(plass_cov)), 
                            Position = 1:nrow(plass_cov), 
                            Coverage=as.numeric(unlist(plass_cov)), 
                            Bin = rep(hubin, nrow(plass_cov))) 
  }
  bin_cov <- sapply(coverage(markers2[markers2$origin == "bin"]), as.numeric)
  if(typeof(nrow(bin_cov)) == "NULL"){
    covdf <- rbind(covdf, plass_cov)
    next
  } else {
    bin_cov <- data.frame(pfam=rep("gyrB", nrow(bin_cov)), 
                          Origin=rep("bin", nrow(bin_cov)), 
                          Position = 1:nrow(bin_cov), 
                          Coverage=as.numeric(unlist(bin_cov)), 
                          Bin = rep(hubin, nrow(bin_cov)))
    covdf_gyrB <- rbind(covdf_gyrB, plass_cov, bin_cov)
  }
}

covdf2_gyrB <- covdf_gyrB %>% filter(Coverage != 0)

ggplot(covdf_gyrB %>% filter(Coverage != 0), aes(Position, Coverage, fill=Origin)) + 
  geom_bar(stat="identity", position="identity") + facet_wrap(~ Bin) +
  theme_minimal() +
  xlim(c(1, 293)) + 
  scale_y_sqrt() + 
  ggtitle("Per Amino Acid Coverage of PFAM gyrB") +
  ylab("Square Root of Coverage Depth")
