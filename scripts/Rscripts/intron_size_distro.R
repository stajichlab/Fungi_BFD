#!/usr/bin/env Rscript

library(tidyverse)
library(duckdb)
library(purrr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(microshades)
library(forcats)
library(scales)
library(colorspace)
library(ggrepel)
library(purrr)


plot_intronlen_distro <- function(dt,max_intron=1000,min_intron=40) {
  dt$FAMILY = factor(dt$FAMILY)
  nb.cols <- length(unique(dt$FAMILY))
  mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
  
  p = ggplot(dt %>% filter(intron_length < max_intron & intron_length > min_intron),  
             aes(x = intron_length, fill = FAMILY, color=FAMILY)) + 
    scale_colour_manual(values = mycolors) +
    scale_fill_manual(values = mycolors) +
    geom_density(alpha = 0.6) + 
    #geom_histogram(bins=100,alpha = 0.6) + 
    theme_cowplot(12) 
  p
}

plot_intronlen_class_distro <- function(dt,max_intron=1000,min_intron=40) {
  dt$CLASS = factor(dt$CLASS)
  nb.cols <- length(unique(dt$CLASS))
  mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
  
  p = ggplot(dt %>% filter(intron_length < max_intron & intron_length > min_intron),  
             aes(x = intron_length, fill = CLASS, color=CLASS)) + 
    scale_colour_manual(values = mycolors) +
    scale_fill_manual(values = mycolors) +
    # (..count..)/sum(..count..))
    #after_stat(count / ave(count, , FUN = sum))
    geom_density(binwidth=2,alpha = 0.6) + theme_cowplot(12) +
    scale_y_continuous(labels = scales::percent)
  p
}

plot_intronlen_order_distro <- function(dt,max_intron=1000,min_intron=40) {
  dt$ORDER = factor(dt$ORDER)
  nb.cols <- length(unique(dt$ORDER))
  mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
  
  p = ggplot(dt %>% filter(intron_length < max_intron & intron_length > min_intron),  
             aes(x = intron_length, fill = ORDER, color=ORDER)) + 
    scale_colour_manual(values = mycolors) +
    scale_fill_manual(values = mycolors) +
    geom_histogram(bins=100,alpha = 0.6) + theme_cowplot(12) 
  p
}

plot_intronlen_distro_facet <- function(dt,max_intron=500,min_intron=40) {
  dt$FAMILY = factor(dt$FAMILY)
  #  nb.cols <- length(unique(dt$ORDER))
  #  mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
  
  p = ggplot(dt %>% filter(intron_length < max_intron & intron_length > min_intron),  
             aes(x = intron_length, fill = FAMILY, color=FAMILY)) + 
    geom_density(alpha = 0.4) + 
    facet_wrap(dt$SUBPHYLUM) + 
    theme_cowplot(12)
  p
}

# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir="intronDB/introns.duckdb", read_only = TRUE)

speciessql = "SELECT s.LOCUSTAG, PHYLUM, SUBPHYLUM, FAMILY, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN, TOTAL_LENGTH, N50
FROM species as s, asm_stats as stats
WHERE s.LOCUSTAG = stats.LOCUSTAG"

speciesinfo <- dbGetQuery(con, speciessql)

# intron length calculated
intronlensql="
SELECT gi.LOCUSTAG, intron_length, PHYLUM, SUBPHYLUM, CLASS, FAMILY, s.ORDER, GENUS, s.SPECIES, s.STRAIN
FROM species s,
 (SELECT substring(transcript_id,1,8) as LOCUSTAG, 
         abs(gene_introns.end - gene_introns.start) as intron_length
  FROM gene_introns) as gi
WHERE gi.LOCUSTAG = s.LOCUSTAG"

intronlens <- dbGetQuery(con, intronlensql)
head(intronlens)

p2 <- plot_intronlen_distro_facet(intronlens,200,30)


Kickxello = intronlens %>% filter(SUBPHYLUM=="Kickxellomycotina")
plot_k <- plot_intronlen_distro(Kickxello,500)
plot_k

Sacch = intronlens %>% filter(SUBPHYLUM=="Saccharomycotina")
plot_s<- plot_intronlen_distro(Sacch,500)
plot_s

plot_s<- plot_intronlen_distro(Sacch,200)
plot_s

plot_s<- plot_intronlen_distro(Sacch,100)
plot_s

Agar = intronlens %>% filter(SUBPHYLUM=="Agaricomycotina")
plot_a<- plot_intronlen_order_distro(Agar,100) + 
  ggtitle("Agaricomycotina intron length distribution by Order")
plot_a

plot_a<- plot_intronlen_distro(Agar,100) + 
  ggtitle("Agaricomycotina intron length distribution by Order")
plot_a

Pez = intronlens %>% filter(SUBPHYLUM=="Pezizomycotina")
plot_p<- plot_intronlen_order_distro(Pez,200) + 
  ggtitle("Pezizomycotina intron length distribution by Order")
plot_p

plot_p<- plot_intronlen_distro(Pez,200) + 
  ggtitle("Pezizomycotina intron length distribution by Order")
plot_p

Muc = intronlens %>% filter(SUBPHYLUM=="Mucoromycotina")
plot_m<- plot_intronlen_distro(Muc,150) + ggtitle("Mucoromycotina intron length distribution by Order")
plot_m

AMF = intronlens %>% filter(SUBPHYLUM=="Glomeromycotina")
plot_amf<- plot_intronlen_distro(AMF,150) + ggtitle("Glomeromycotina intron length distribution by Order")
plot_amf

Taph = intronlens %>% filter(SUBPHYLUM=="Taphrinomycotina")
plot_tap<- plot_intronlen_distro(Taph,70,25) + ggtitle("Taphrinomycotina intron length distribution by Order")
plot_tap

Chytrid = intronlens %>% filter(SUBPHYLUM=="Chytridiomycotina")
plot_chytrid <- plot_intronlen_distro(Chytrid,150,25) + ggtitle("Chytridiomycotina intron length distribution by Order")
plot_chytrid

Blasto = intronlens %>% filter(SUBPHYLUM=="Blastocladiomycotina")
plot_Bchytrid <- plot_intronlen_class_distro(Blasto,150,25) + 
  ggtitle("Blastocladiomycotina intron length distribution by Class")
plot_Bchytrid


Mic = intronlens %>% filter(SUBPHYLUM=="Microsporidiomycotina")
plot_mic<- plot_intronlen_distro(Mic,500,1) + ggtitle("Microsporidiomycotina intron length distribution by Order")
plot_mic

dbDisconnect(con, shutdown = TRUE)


