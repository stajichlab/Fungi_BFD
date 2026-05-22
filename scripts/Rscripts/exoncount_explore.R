#!/usr/bin/env Rscript

library(tidyverse)
library(duckdb)
library(purrr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

plotdir = "plots/exon_explore"
dir.create(file.path(plotdir),showWarnings = FALSE)
pdf(file.path(plotdir,"exon_explore.pdf"))
# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir="functionalDB/function.duckdb", read_only = TRUE)

exoncount_sql = "SELECT species.PHYLUM, species.SUBPHYLUM, species.GENUS, species.SPECIES, asm_stats.TOTAL_LENGTH, 
                 asm_stats.GC_PERCENT, ex.*, genec.genecount
FROM 
(select ge.LOCUSTAG, ge.transcript_id, count(*) as exoncount FROM gene_exons ge, gene_transcripts gt, gene_info gi 
WHERE 
 gi.gene_id = gt.gene_id AND gt.transcript_id = ge.transcript_id AND 
 gi.gene_type = 'protein_coding' GROUP BY ge.LOCUSTAG, ge.transcript_id) as ex,
 (select gt.LOCUSTAG, count(*) as genecount FROM gene_transcripts gt GROUP BY LOCUSTAG) as genec,
 species, asm_stats
WHERE ex.LOCUSTAG = species.LOCUSTAG AND ex.LOCUSTAG = asm_stats.LOCUSTAG and genec.LOCUSTAG = ex.LOCUSTAG"

exons_sumstat <- dbGetQuery(con, exoncount_sql) 
exons_sumstat
exons_sumstatfilter <- exons_sumstat %>% filter(! PHYLUM %in% c("Sanchytriomycota","Cryptomycota","NA") & ! is.na(PHYLUM) ) %>% 
  filter(exoncount <= 10) %>% mutate(exoncount = as_factor(exoncount))


mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(exons_sumstatfilter$SUBPHYLUM)))
myphylumcolors <- colorRampPalette(brewer.pal(7, "Set1"))(length(unique(exons_sumstatfilter$PHYLUM)))

p <- ggplot(exons_sumstatfilter,
            aes(x=exoncount, y= after_stat(count), color=SUBPHYLUM, fill=SUBPHYLUM)) +
  #geom_histogram(alpha=0.5, stat="count",binwidth=1) +
  geom_bar(alpha=0.5) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Density") +
  xlab("Transcript Exon Count") + ggtitle("Histogram of Exon counts") +
  theme_cowplot(12) + theme(legend.position="bottom") + facet_wrap(~PHYLUM,scales = 'free_y')
  
p

ggsave(file.path(plotdir,"exon_count_histogram.pdf"),p,width=10,height=8)

exoncountconserved_sql = "SELECT species.PHYLUM, species.SUBPHYLUM, species.GENUS, species.SPECIES, ex.*, og.orthogroup, ogc.og_size
FROM 
(select ge.LOCUSTAG, ge.transcript_id, count(*) as exoncount FROM gene_exons ge, gene_transcripts gt, gene_info gi 
WHERE 
 gi.gene_id = gt.gene_id AND gt.transcript_id = ge.transcript_id AND 
 gi.gene_type = 'protein_coding' GROUP BY ge.LOCUSTAG, ge.transcript_id
) as ex, mmseqs_orthogroup_clusters og,
(SELECT orthogroup, count(*) as og_size FROM mmseqs_orthogroup_clusters GROUP BY orthogroup) as ogc,
 species
WHERE ex.LOCUSTAG = species.LOCUSTAG AND ex.transcript_id = og.transcript_id AND 
ogc.orthogroup = og.orthogroup AND ogc.og_size > 1"

exonsconserved_sumstat <- dbGetQuery(con, exoncountconserved_sql) 
exonsconserved_sumstat

exonsconserved_sumstat <- exonsconserved_sumstat %>% 
  filter(! PHYLUM %in% c("Sanchytriomycota","Cryptomycota","NA") & ! is.na(PHYLUM) ) %>% 
  filter(exoncount <= 10) %>% 
  mutate(exoncount = as.factor(exoncount))

mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(exonsconserved_sumstat$SUBPHYLUM)))
myphylumcolors <- colorRampPalette(brewer.pal(7, "Set1"))(length(unique(exonsconserved_sumstat$PHYLUM)))

p <- ggplot(exonsconserved_sumstat,
            aes(x=exoncount, y= after_stat(count), color=SUBPHYLUM, fill=SUBPHYLUM)) +
  #geom_histogram(alpha=0.5, stat="count",binwidth=1) +
  geom_bar(alpha=0.5) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Density") +
  xlab("Transcript Exon Count") + ggtitle("Histogram of Exon counts With Orthologs") +
  theme_cowplot(12) + theme(legend.position="bottom") + facet_wrap(~PHYLUM,scales = 'free_y')

p

ggsave(file.path(plotdir,"exon_count_OnlyConserved_histogram.pdf"),p,width=10,height=8)


#
category_exoncountconserved_sql = "
SELECT species.PHYLUM, species.SUBPHYLUM, species.GENUS, species.SPECIES, 
categoryexonct.*, gc.genecount, categoryexonct.count_exon_category / gc.genecount as exoncategory_abundance,
FROM
(SELECT ex.LOCUSTAG, ex.exoncount, COUNT(*) as count_exon_category
FROM 
(select ge.LOCUSTAG, ge.transcript_id, count(*) as exoncount FROM gene_exons ge, gene_transcripts gt, gene_info gi 
WHERE 
 gi.gene_id = gt.gene_id AND gt.transcript_id = ge.transcript_id AND 
 gi.gene_type = 'protein_coding' GROUP BY ge.LOCUSTAG, ge.transcript_id
) as ex, mmseqs_orthogroup_clusters og,
(SELECT orthogroup, count(*) as og_size FROM mmseqs_orthogroup_clusters GROUP BY orthogroup) as ogc
WHERE ex.transcript_id = og.transcript_id AND 
ogc.orthogroup = og.orthogroup AND ogc.og_size > 1
GROUP BY ex.LOCUSTAG, ex.exoncount) as categoryexonct,
(select gt.LOCUSTAG, count(*) as genecount FROM gene_transcripts gt GROUP BY LOCUSTAG) as gc,
species
WHERE categoryexonct.LOCUSTAG = species.LOCUSTAG AND gc.LOCUSTAG=categoryexonct.LOCUSTAG"

category_exoncountconserved_sumstat <- dbGetQuery(con, category_exoncountconserved_sql) 
category_exoncountconserved_sumstat

category_exoncountconserved_sumstatfilter <- category_exoncountconserved_sumstat %>% 
  filter(! PHYLUM %in% c("Sanchytriomycota","Cryptomycota","NA") & ! is.na(PHYLUM) ) %>% 
  filter(exoncount <= 10) %>% 
  mutate(exoncount = as.factor(exoncount))


mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(category_exoncountconserved_sumstatfilter$SUBPHYLUM)))
myphylumcolors <- colorRampPalette(brewer.pal(7, "Set1"))(length(unique(category_exoncountconserved_sumstatfilter$PHYLUM)))

p <- ggplot(category_exoncountconserved_sumstatfilter,
            aes(x=exoncount, y=exoncategory_abundance, color=SUBPHYLUM, fill=SUBPHYLUM)) +
  geom_point(alpha=0.5) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Density") +
  xlab("Transcript Exon Count") + ggtitle("Frequency of transcripts with exon counts with Orthologs") +
  theme_cowplot(12) + theme(legend.position="bottom") + facet_wrap(~PHYLUM,scales = 'free_y')

p
ggsave(file.path(plotdir,"transcript_exoncount_OnlyConserved_histogram.pdf"),p,width=10,height=8)

dbDisconnect(con, shutdown = TRUE)


