mh08blast <- read.table("outputs/hu-croissants/assembly/blast/hu-genome08-blastp.tab")

# Grab highest bitscore match 
mh08blast<- mh08blast[order(mh08blast[, 1], mh08blast[ , 12]), ]

# Select the last row by id
mh08blast_bitscore <- mh08blast[!duplicated(mh08blast[ , 1], fromLast=TRUE), ]
colnames(mh08blast_bitscore) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                   "qstart", "qend", "sstart", "send", "evalue", "bitscore")


mh41blast <- read.table("outputs/hu-croissants/assembly/blast/hu-genome41-blastp.tab")

# Grab highest bitscore match 
mh41blast<- mh41blast[order(mh41blast[, 1], mh41blast[ , 12]), ]

# Select the last row by id
mh41blast_bitscore <- mh41blast[!duplicated(mh41blast[ , 1], fromLast=TRUE), ]
colnames(mh41blast_bitscore) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                                  "qstart", "qend", "sstart", "send", "evalue", "bitscore")
