#!/usr/bin/env Rscript

library(tidyverse)
library(duckdb)
library(purrr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

pdf("plots/genome_size_explore.pdf")

# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir="intronDB/introns.duckdb", read_only = TRUE)

genecountsql ="
SELECT sp.LOCUSTAG, sp.GENUS, sp.SPECIES, pg.gene_count
FROM (SELECT LOCUSTAG, COUNT(*) as gene_count FROM 
(SELECT substring(transcript_id,1,8) as LOCUSTAG FROM gene_proteins p)
GROUP BY LOCUSTAG) as pg,
species as sp
WHERE sp.LOCUSTAG = pg.LOCUSTAG"

genecount <- dbGetQuery(con, genecountsql)

target_genera = c('Aspergillus','Penicillium','Cryptococcus','Coccidioides','Histoplasma','Sporothrix','Paracoccidioides')
targets = tibble(genus = target_genera)

temptablesql = "DROP TABLE IF EXISTS target_genera; CREATE TEMPORARY TABLE target_genera (genus VARCHAR);"
dbExecute(con, temptablesql)
dbAppendTable(con,"target_genera",targets)

temptable = "SELECT * from target_genera"
q = dbGetQuery(con, temptable)
q

asmstats_sql = "SELECT sp.SUBPHYLUM, sp.GENUS, asm_stats.* 
FROM asm_stats, species sp, target_genera
WHERE 
sp.LOCUSTAG = asm_stats.LOCUSTAG AND sp.GENUS = target_genera.genus"
asm_stats <- dbGetQuery(con, asmstats_sql)

# barplot for these GENERA

sumsizes <- asm_stats %>% group_by(GENUS) %>% summarise(
                                g_mean = mean(TOTAL_LENGTH / 1000000),
                                g_N = n(),
                                g_sd = sd(TOTAL_LENGTH/1000000),
                                g_se = g_sd / sqrt(g_N))


nb.cols <- 8
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

p <- ggplot() + geom_bar(data=sumsizes,
                           aes(y=g_mean,
                               x=GENUS,
                               fill=GENUS,
                               color=GENUS),
                           stat="identity", width=0.5) +
    geom_errorbar(data = sumsizes,
                  aes(y=g_mean,
                      x=GENUS,
                      ymin = g_mean - g_se,
                      ymax = g_mean + g_se,
                      color = GENUS),
                  stat="identity", width=0.5) +
  geom_point(data=asm_stats, aes(y=TOTAL_LENGTH/1000000, x=GENUS, fill=GENUS),size=2) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  xlab("Genus") +
  ylab("Genome size (Mb)") +
  theme_cowplot(12) 
  
p

ggsave("plots/genome_size_Brazilset.pdf",p,width=10,height=10)



asmstats_sql = "SELECT sp.PHYLUM,sp.SUBPHYLUM, sp.GENUS, asm_stats.* 
FROM asm_stats, species sp
WHERE 
sp.LOCUSTAG = asm_stats.LOCUSTAG"
asm_stats <- dbGetQuery(con, asmstats_sql)

# barplot for these GENERA

sumsizes <- asm_stats %>% group_by(SUBPHYLUM,PHYLUM) %>% summarise(
  g_mean = mean(TOTAL_LENGTH / 1000000),
  g_N = n(),
  g_sd = sd(TOTAL_LENGTH/1000000),
  g_se = g_sd / sqrt(g_N)) %>% mutate(SUBPHYLUM = fct_reorder(SUBPHYLUM, PHYLUM))


nb.cols <- 23
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

p <- ggplot() + geom_bar(data=sumsizes,
                         aes(y=g_mean,
                             x=SUBPHYLUM,
                             fill=SUBPHYLUM,
                             color=SUBPHYLUM),
                         stat="identity", width=0.5) +
  geom_errorbar(data = sumsizes,
                aes(y=g_mean,
                    x=SUBPHYLUM,
                    ymin = g_mean - g_se,
                    ymax = g_mean + g_se,
                    color = SUBPHYLUM),
                stat="identity", width=0.5) +
  geom_point(data=asm_stats, aes(y=TOTAL_LENGTH/1000000, x=SUBPHYLUM, fill=SUBPHYLUM),size=2) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  xlab("Subphylum") +
  ylab("Genome size (Mb)") +
  theme_cowplot(12) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14)) + scale_y_log10() +
  theme(legend.position="none")

p

ggsave("plots/genome_size_barplot.pdf",p,width=10,height=10)

asmstats_sql = "SELECT sp.PHYLUM,sp.SUBPHYLUM, sp.GENUS, asm_stats.*, gene_count
FROM asm_stats, species sp,
(SELECT LOCUSTAG, count(*) as gene_count
FROM gene_proteins 
GROUP BY LOCUSTAG) as gp
WHERE gp.LOCUSTAG = sp.LOCUSTAG AND sp.LOCUSTAG = asm_stats.LOCUSTAG"
asm_stats <- dbGetQuery(con, asmstats_sql)

nb.cols <- 23
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

p <- ggplot(asm_stats) +
  geom_point(aes(x=TOTAL_LENGTH/1000000, y=gene_count, color=SUBPHYLUM, fill=SUBPHYLUM),size=2) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Gene count") +
  xlab("Genome size (Mb)") + ggtitle("Gene Count vs Genome Size") +
  theme_cowplot(12) + theme(legend.position="bottom")
p
ggsave("plots/genome_size_by_subphylum.pdf",p,width=15,height=10)

dbDisconnect(con, shutdown = TRUE)


