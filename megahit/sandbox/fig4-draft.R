library(tidyverse)
setwd("~/github/hu-snake/")
blast <- read.csv("sandbox/info_blast_small.csv", stringsAsFactors = F)
blast$sb <- c("SB2", "SB2", "SB2", "SB2", "SB2", "SB1", "SB1", "SB1","SB1", "SB1", "SB1","SB1", "SB1", "SB1")
blast$color <- c("black", "black", "black", "black", "black", "red", "red", "red","red", "red", "red","red", "red", "red")
blast$hu_taxonomy <- trimws(blast$hu_taxonomy)

blast_sb1 <- blast %>% filter(sb =="SB1")
blast_sb2 <- blast %>% filter(sb =="SB2")
sb1 <- ggplot(blast_sb1) + 
  geom_point(aes(x = AA.PIdent, y = as.numeric(row.names(blast_sb1)), col = Marker.Gene), size = 3) +
  theme_bw() +
  facet_wrap(~sb) +
  scale_y_continuous(labels = as.character(blast_sb1$hu_taxonomy), breaks = 1:9, name = "Query",
                     sec.axis = dup_axis(labels = as.character(blast_sb1$AA.Species), 
                                         breaks = 1:9, name = "Matches"))
sb2 <- ggplot(blast_sb2) + 
  geom_point(aes(x = AA.PIdent, y = as.numeric(row.names(blast_sb2)), col = Marker.Gene), size = 3) +
  theme_bw() +
  facet_wrap(~sb) +


gridExtra::grid.arrange(sb1, sb2, nrow=2)

blast$AA.Species2 <- ifelse(blast$AA.Species == blast$hu_taxonomy, "", blast$AA.Species)

blast <- blast[order(blast$AA.PIdent), ]

plt <- ggplot(blast) + 
        geom_point(aes(x = AA.PIdent, y = 1:nrow(blast), col = Marker.Gene, size = AA.Length)) +
        theme_bw() +
        scale_y_continuous(labels = as.character(paste0(blast$hu_taxonomy, " (", blast$sb, ")")), breaks = 1:nrow(blast), name = "Query",
                            sec.axis = dup_axis(labels = as.character(blast$AA.Species2), 
                                         breaks = 1:nrow(blast), name = "Matches")) +
        scale_x_continuous(limits = c(89, 100))
        #theme(axis.text.y = element_text(colour = blast$color))
plt

