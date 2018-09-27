library(GenomicRanges)
library(dplyr)

source(snakemake@input[['rhmmer']])

dom <- read_domtblout(snakemake@input[['dom']])

# read in info about gyrA seqs that matched pfam domain
dom <- makeGRangesFromDataFrame(dom,
                                keep.extra.columns=T,
                                ignore.strand=T,
                                seqinfo=NULL,
                                seqnames.field=c("domain_name"),
                                start.field="hmm_from",
                                end.field="hmm_to",
                                #strand.field="strand",
                                starts.in.df.are.0based=FALSE)

# define 200 bp interval to calc pid with based on max coverage
cov <- coverage(dom)
mid <- which.max(cov)
lower <- mid - 100
upper <- mid + 100

# define reads to keep for pid calc ------------------------------------

# define domain_name using dom object
# start and end based on cov profile

# make a sliding window of 200bp that progresses by 10 bp each step. Calc overlaps for each of these. Select the one that is longest at the end. 
ranges <- seq(from = 1, to = (max(dom$domain_len) - 200), by = 10)
all <- list()

for(start in ranges){
  end = start + 200
  want <- makeGRangesFromDataFrame(data.frame(domain_name = as.character(seqnames(dom))[1], 
                                              start = start, end = end), 
                                   keep.extra.columns=T,
                                   ignore.strand=T,
                                   seqinfo=NULL,
                                   seqnames.field=c("domain_name"),
                                   start.field="start",
                                   end.field="end",
                                   #strand.field="strand",
                                   starts.in.df.are.0based=FALSE)
  ov <- subsetByOverlaps(dom, want, minoverlap = 125)
  all[[length(all)+1]] <- ov
}

# find the granges object with the most records
most <- which.max(sapply(all, length))

print(paste0("200 aa region with most overlapping seqs starts at residue ", ranges[most]))

# write sequence names to a file
write.table(all[[most]]$query_name, file = snakemake@output[['keep']], quote = F, row.names = F, col.names = F)
