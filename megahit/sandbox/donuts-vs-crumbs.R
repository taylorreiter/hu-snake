setwd("github/hu-snake/inputs/hu-s1_k31_r1_search_oh0_may20/")
library(Biostrings)
library(tidyverse)
library(ggthemes)

# Compare Unitigs (number & nucleotides) ----------------------------------

# Read in data
    donut_files <- list.files(path = ".", pattern = ".contigs.fa.gz.donut.fa$")
    crumb_files <- list.files(path = ".", pattern = ".contigs.fa.gz.crumbs.fa$")
    
    donuts <- lapply(donut_files, readDNAStringSet)
    crumbs <- lapply(crumb_files, readDNAStringSet)

# Accuracy check. Should print TRUE.
# Is the number of unitigs in the crumbs equal to the number of unitigs from the crumbs contained in the donuts?
    for(i in 1:length(crumbs)){
      crumb <- crumbs[[i]]
      donut <- donuts[[i]]
      print(length(crumb[crumb %in% donut]) == length(crumb))
    }


# how many unitigs are subtracted from the donuts to produce the crumbs?
    crumb_uni_total <- rep(NA, length(crumbs))
    donut_uni_total <- rep(NA, length(crumbs))
    donut_uni_no_crumbs <- rep(NA, length(crumbs))
    for(i in 1:length(crumbs)){
      crumb <- crumbs[[i]]
      donut <- donuts[[i]]
      crumb_uni_total[i] <- length(crumb)
      donut_uni_total[i] <- length(donut)
      donut_uni_no_crumbs[i] <- length(donut) - length(crumb)
    }
    
    # plot of total donut unitigs, colored by amount subtracted to make crumbs
    unis <- as.data.frame(donut_uni_no_crumbs)
    unis$crumb_uni_total  <- crumb_uni_total 
    # label df rows
    unis$bin <- names
    unis <- gather(unis, key = uni_type, value = number, -bin)
    ggplot(unis, aes(x = bin, y = number, fill = uni_type)) + 
      theme_pander() +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle("Number of Unitigs in Donuts") +
      xlab("Query Bin") +
      ylab("Number of unitigs") +
      labs(fill = "Donut Composition") +
      scale_fill_hue(labels = c("Crumb & Donut", "Donut Only")) 

# how many more nucleotides are there in the donut unitigs than in the crumbs unitigs?
    crumb_nucs_total <- rep(NA, length(crumbs))
    donut_nucs_total <- rep(NA, length(crumbs))
    donut_nucs_no_crumbs <- rep(NA, length(crumbs))
    for(i in 1:length(crumbs)){
      crumb <- crumbs[[i]]
      donut <- donuts[[i]]
      donut_nucs_total[i] <- sum(nchar(donut))
      crumb_nucs_total[i] <- sum(nchar(crumb))
      donut_nucs_no_crumbs[i] <- sum(nchar(donut)) - sum(nchar(donut[donut %in% crumb]))
    }
    
    # plot of total donut nucleotides, colored by amount subtracted to make crumbs
    nucs <- as.data.frame(donut_nucs_no_crumbs)
    nucs$donut_nucs_total <- donut_nucs_total
    nucs$crumb_nucs_total <- crumb_nucs_total
    # label df rows
    names <- gsub(".fa.cdbg_ids.contigs.fa.gz.donut.fa", "", donut_files)
    nucs$bin <- names
    nucs <- gather(nucs, key = nuc_type, value = number, -bin)
    ggplot(nucs %>% filter(nuc_type != "donut_nucs_total"), aes(x = bin, y = number, fill = nuc_type)) + 
      theme_pander() +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle("Number of Nucleotides in Donuts") +
      xlab("Query Bin") +
      ylab("Number of nucleotides (bp)") +
      labs(fill = "Donut Composition") +
      scale_fill_hue(labels = c("Crumb & Donut", "Donut Only")) 

# What is the distribution of sizes of the unitigs that are in the donuts but not in the crumbs??
    dist <- list()
    for(i in 1:length(crumbs)){
      crumb <- crumbs[[i]]
      donut <- donuts[[i]]
      dist[[i]] <- summary(nchar(donut[!(donut %in% crumb)]))
    }

    
donuts_all <- do.call("c", donuts)
donuts_char <- as.data.frame(nchar(donuts_all))
ggplot(donuts_char, aes(x = nchar(donuts_all))) + geom_histogram() +
  theme_pander() +
  scale_y_log10() +
  ylab("log10(count)") +
  xlab("length in nucleotides") +
  ggtitle("Length of unitigs in donuts")

crumbs_all <- do.call("c", crumbs)
crumbs_char <- as.data.frame(nchar(crumbs_all))
ggplot(crumbs_char, aes(x = nchar(crumbs_all))) + geom_histogram() +
  theme_pander() +
  scale_y_log10() +
  ylab("log10(count)") +
  xlab("length in nucleotides") +
  ggtitle("Length of unitigs in crumbs")

# Compare Reads (number & nucleotides) ------------------------------------
    


    # Read in data
    donut_reads_files <- list.files(path = ".", pattern = ".reads.fa.gz.donut.fa$")
    crumb_reads_files <- list.files(path = ".", pattern = ".reads.fa.gz.crumbs.fa$")
    
    donuts_reads <- lapply(donut_reads_files, readDNAStringSet)
    crumbs_reads <- lapply(crumb_reads_files, readDNAStringSet)
    
    # Accuracy check. Should print TRUE.
    # Is the number of unitigs in the crumbs equal to the number of unitigs from the crumbs contained in the donuts?
    for(i in 1:length(crumbs_reads)){
      crumb <- crumbs_reads[[i]]
      donut <- donuts_reads[[i]]
      print(length(crumb[crumb %in% donut]) == length(crumb))
    }
    
    
    # how many unitigs are subtracted from the donuts to produce the crumbs?
    crumb_reads_total <- rep(NA, length(crumbs_reads))
    donut_reads_total <- rep(NA, length(crumbs_reads))
    donut_reads_no_crumbs <- rep(NA, length(crumbs_reads))
    for(i in 1:length(crumbs_reads)){
      crumb <- crumbs_reads[[i]]
      donut <- donuts_reads[[i]]
      crumb_reads_total[i] <- length(crumb)
      donut_reads_total[i] <- length(donut)
      donut_reads_no_crumbs[i] <- length(donut) - length(crumb)
    }
    
    # plot of total donut unitigs, colored by amount subtracted to make crumbs
    reads <- as.data.frame(donut_reads_no_crumbs)
    reads$crumb_reads_total  <- crumb_reads_total 
    # label df rows
    names <- gsub(".fa.cdbg_ids.reads.fa.gz.donut.fa", "", donut_reads_files)
    reads$bin <- names
    reads <- gather(reads, key = read_type, value = number, -bin)
    ggplot(reads, aes(x = bin, y = number, fill = read_type)) + 
      theme_pander() +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle("Number of Reads in Donuts") +
      xlab("Query Bin") +
      ylab("Number of Reads") +
      labs(fill = "Donut Composition") +
      scale_fill_hue(labels = c("Crumb & Donut", "Donut Only"))
    