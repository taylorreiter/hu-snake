best <- read.delim("explore/prots/best-archaea-gyra-blastp.tab", sep = "\t", header = F, 
                            col.names = c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 
                                          'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'))
best$bin <- gsub("(_SRR1976948.[0-9_]{2,12})", "", best$qseqid)
