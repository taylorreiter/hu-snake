setwd("~/github/hu-snake/sandbox/PLASS/explore/markers/hmmer/")

## In Bash

## check if this is correct
# conda install -c bioconda hmmer3...

## download hmm
# wget -O PF00521_full_gyra.sto https://pfam.xfam.org/family/PF00521/alignment/full
#
## make hmm
# hmmbuild PF00521_full_gyra.hmm PF00521_full_gyra.txt
## locate the gyra sequences in the plass all_sb1 output. Use a bitscore cutoff of 100.
# hmmscan -T 100 -o plass-gyra-hmmscanT100.out --tblout plass-gyra-seq-hmmscanT100.out --domtblout plass-gyra-dom-hmmscanT100.out --pfamtblout plass-gyra-pfam-hmmscanT100.out PF00521_full_gyra.hmm all_sb1.nostop.c100.fas 
# 
## grab names of sequences that matched (there are 2170)
# grep -v "^#" plass-gyra-seq-hmmscanT100.out | cut -d " " -f15 > plass-gyra-seq-hmmscanT100-NAMES.out
# 
## grab the fasta sequences
# while read inline
# do 
# samtools faidx plass-gyra-hmmscan.faa $inline >> plass-gyra-hmmscanT100.faa
# done < plass-gyra-seq-hmmscanT100-NAMES.out
#
## align all of the fasta sequences, using the hmm model
# hmmalign -o PF00521_plass_gyra.sto PF00521_full_gyra.hmm plass-gyra-hmmscanT100.faa 

library(tidyr)
library(dplyr)
library(plotly)
library(GenomicRanges)
library(rhmmer)
library(ggplot2)
library(ggbio)

# coverage plot to see what aas to limit to -------------------------------

dom <- read_domtblout('plass-gyra-dom-hmmscanT100.out')

dom <- dom %>%
  mutate(bin = gsub("(_SRR1976948.[0-9_]{2,12})", "", query_name)) %>% 
  mutate(gene = rep("gyrA", nrow(dom))) 

# plot
dom <- makeGRangesFromDataFrame(dom,
                                keep.extra.columns=T,
                                ignore.strand=T,
                                seqinfo=NULL,
                                seqnames.field=c("domain_name"),
                                start.field="hmm_from",
                                end.field="hmm_to",
                                #strand.field="strand",
                                starts.in.df.are.0based=FALSE)

ggplot(dom) + 
  stat_coverage(geom = "bar") +
  theme_minimal() +
  scale_y_sqrt() 

# limit to 200-400 --------------------------------------------------------

want <- makeGRangesFromDataFrame(data.frame(domain_name = "PF00521_full_gyra", start = 200, end = 400), 
                                 keep.extra.columns=T,
                                 ignore.strand=T,
                                 seqinfo=NULL,
                                 seqnames.field=c("domain_name"),
                                 start.field="start",
                                 end.field="end",
                                 #strand.field="strand",
                                 starts.in.df.are.0based=FALSE)


ov <- subsetByOverlaps(dom, want, minoverlap = 125)
table(ov$bin) # check that all bins are well represented
#write.table(ov$query_name, file = "plass-gyra-hmmscanT100-gyra200-400-NAMES.txt", 
#            quote = F, row.names = F, col.names = F)

## In bash
# while read inline
# do
# samtools faidx plass-gyra-hmmscanT100.faa $inline >> plass-gyra-hmmscanT100-gyra200-400.faa
# done < plass-gyra-hmmscanT100-gyra200-400-NAMES.txt 
# 
## align these new sequences using hmm
# hmmalign -o PF00521_plass_gyra200-400.sto PF00521_full_gyra.hmm plass-gyra-hmmscanT100-gyra200-400.faa 
# esl-alipid PF00521_plass_gyra200-400.sto > PF00521_plass_gyra200-400.sto.pid

# plot --------------------------------------------------------------------

#pid <- read.table(file = "hmmalign/PF00521_plass_gyra200-400.sto.pid", skip = 1, header = F)
## clustal alignmnet of these seqs was better. 
## clustal was done online
## In the future:
#  conda install -c bioconda clustalw 

## Then in bash:
# esl-alimanip --informat clustal --outformat stockholm CLUSTALAplass-gyra-hmmscanT100-gyra200-400.clu_num # or something similar...
# esl-alipid CLUSTALAplass-gyra-hmmscanT100-gyra200-400.sto > CLUSTALAplass-gyra-hmmscanT100-gyra200-400.sto.pid
pid <- read.table(file = "hmmalign/CLUSTALAplass-gyra-hmmscanT100-gyra200-400.sto.pid", skip = 1, header = F)
colnames(pid) <- c('seqname1', 'seqname2', 'pid', 'nid', 'denomid', 'pmatch', 'nmatch', 'denommatch')
hist(pid$pid)
pid$alipid <- pid$nid / pid$nmatch # calculate percent identity only over length of sequence where two sequences overlapped
pid$alipid <- 1-pid$alipid # change into distance, instead of similarity
pid2 <- pid[ , c('seqname1', 'seqname2', 'alipid')]
# make sure all combinations are represented
dummy1 <- data.frame(seqname1 = pid2$seqname2, seqname2 = pid2$seqname1, alipid = pid2$alipid)
dummy2 <- data.frame(seqname1 = unique(pid2$seqname1), seqname2 = unique(pid2$seqname1), alipid = rep(0, length(unique(pid2$seqname1)))) # add in self comp for diag of matrix
pid2 <- rbind(pid2, dummy1, dummy2)

pidw <- spread(pid2, key = seqname2, value = alipid, drop = F)
row.names(pidw) <- pidw[ , 'seqname1']
pidw <- pidw[ , -1]
pidw <- as.matrix(pidw)
pidw[is.na(pidw)] <- 1
mds <- cmdscale(pidw)

mdsdf <- as.data.frame(mds)
mdsdf$name <- gsub("_.*", "", row.names(mdsdf))

# merge with info to get abbreviated species name
info <- read.csv("../../../inputs/hu_info.csv")
mdsdf <- left_join(mdsdf, info, by = "name")

# Remove anything that occurs less than 7 times
tmp <- mdsdf %>%
        filter(name != "hu-genome19") %>% # 3
        filter(name != "hu-genome23") %>% # 7
        filter(name != "hu-genome37") # 4

# define colors to match Felix
viridian         = "#4b9179ff" # Green with a hint of blue
selective_yellow = "#ffb500ff" # Excavator yellow
rose             = "#ff1164ff"
hot_pink         = "#ff70b0ff" 
steel_blue       = "#3e80b8ff" 
# bone             = "#ded5baff" # Off-white
# white_rock       = "#f1ece0ff" # Lighter off-white
fuscous_gray     = "#5c594fff" # Darkish gray
newblack         = "#0e3452ff" # Almost black, blue-tinted
purple = "#5906a1"
red = "#cc0000"
aqua = "#ace9e7"
mauve = "#b889b6"
orange = "#f68902"
light_brown = "#ba8f80"
bright_green = "#23f210"
grey_yellow = "#6b6b67"

p <- plot_ly(tmp, x = ~V1, y = ~V2,
             # Hover text:
             text = tmp$name,
             color = tmp$name,
             opacity = .75
)
p

ggplot(tmp, aes(x = V1, y = V2, group = abbreviation2)) +
  geom_point(aes(color = abbreviation2), alpha = .75) + theme_minimal() +
  xlab("NMDS1") + ylab("NMDS2") + 
  #ggtitle("MDS Plot of Pairwise Distance of gyrA Percent Identity") +
  scale_color_manual(name = "Bin",
                    values = c(light_brown, selective_yellow, rose, grey_yellow, newblack, 
                                fuscous_gray, steel_blue, mauve, red, aqua, purple, 
                                bright_green, hot_pink, orange, viridian)) +
  theme(legend.text = element_text(face = "italic"))
  
