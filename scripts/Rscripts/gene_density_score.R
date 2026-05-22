#!/usr/bin/env Rscript

library(tidyverse)
library(duckdb)
library(purrr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(paletteer)
library(cowplot)
library(ggfortify)
library(ggpubr)

DBDIR="functionalDB"
DBNAME="function.duckdb"
plotdir = "plots"
dir.create(file.path(plotdir),showWarnings = FALSE)
statsplotdir = file.path(plotdir,"gene_density")
dir.create(statsplotdir, showWarnings = FALSE)
pdf(file.path(plotdir,"gene_density","density_misc.pdf"))
#pdf(file.path(statsplotdir,"misc_plots.pdf"))
# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir=file.path(DBDIR,DBNAME), read_only=TRUE)

densitystat_sql ="
SELECT sp.*, stats.GC_PERCENT, stats.TOTAL_LENGTH, gene_count, mean_gene_length,
 pw_dist.chrom, pw_dist.mean_dist
FROM 
species as sp,
asm_stats as stats,
(SELECT LOCUSTAG, count(*) as gene_count, MEAN(length) as mean_gene_length
FROM gene_proteins 
GROUP BY LOCUSTAG) as gp,
(SELECT species_prefix, chrom, MEAN(distance) as mean_dist
 FROM gene_pairwise_distances as pw_dist, gene_transcripts as gt
 WHERE pw_dist.left_gene = gt.gene_id
GROUP BY species_prefix, chrom) as pw_dist
WHERE 
sp.LOCUSTAG = stats.LOCUSTAG and gp.LOCUSTAG = sp.LOCUSTAG
and pw_dist.species_prefix = sp.LOCUSTAG 
ORDER by sp.LOCUSTAG, chrom"

densitystat_res <- dbGetQuery(con, densitystat_sql)
head(densitystat_res)

hist(subset(densitystat_res$mean_dist,densitystat_res$mean_dist<5000),100)

nb.cols <- length(unique(densitystat_res$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

cdsdensitystat_sql ="
SELECT gcount.chrom_name, ci.length AS chrom_length, 
       gcount.gene_count,
       1000 * (gcount.gene_count / ci.length) as coding_density, 
       sp.*
FROM 
species as sp, chrom_info ci,
(SELECT chrom_info.LOCUSTAG, chrom_info.chrom_name, COUNT(gene_info.*) as gene_count
FROM 
chrom_info, gene_info
WHERE 
gene_info.LOCUSTAG = chrom_info.LOCUSTAG AND
chrom_info.chrom_name = gene_info.chrom 
GROUP BY chrom_info.LOCUSTAG, chrom_info.chrom_name) as gcount
WHERE sp.LOCUSTAG = ci.LOCUSTAG AND 
      gcount.LOCUSTAG = ci.LOCUSTAG AND
      gcount.chrom_name = ci.chrom_name and 
      gene_count > 10
"

cdsdensitystat_res <- dbGetQuery(con, cdsdensitystat_sql)
head(cdsdensitystat_res)

hist(subset(cdsdensitystat_res$coding_density,cdsdensitystat_res$coding_density < 1),100)

nb.cols <- length(unique(cdsdensitystat_res$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

codingdensity_small_p <- ggplot(cdsdensitystat_res %>% filter(chrom_length < 100000)) +
  geom_point(aes(x=chrom_length, y=gene_count,color=SUBPHYLUM, fill=SUBPHYLUM,size=coding_density),alpha=0.5) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Gene Count") +
  xlab("Length (Mb)") + ggtitle("Gene Count vs Chromosome length") +
  theme_cowplot(12) + theme(legend.position="bottom") + facet_wrap(~PHYLUM) + scale_x_log10() 
codingdensity_small_p 

ggsave(file.path(statsplotdir,"genedensity_xyplot_less100k.pdf"),codingdensity_small_p,width=15,height=10)

codingdensity_large_p <- ggplot(cdsdensitystat_res %>% filter(chrom_length >= 100000)) +
  geom_point(aes(x=chrom_length, y=gene_count,color=SUBPHYLUM, fill=SUBPHYLUM,size=coding_density),alpha=0.5) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Gene Count") +
  xlab("Length (Mb)") + ggtitle("Gene Count vs Chromosome length") +
  theme_cowplot(12) + theme(legend.position="bottom") + facet_wrap(~PHYLUM) + scale_x_log10() 
codingdensity_large_p 

ggsave(file.path(statsplotdir,"genedensity_xyplot_more100k.pdf"),codingdensity_large_p,width=15,height=10)

highdensity <- cdsdensitystat_res %>% filter(coding_density >= 0.85) %>% 
  select(c(chrom_name,coding_density,chrom_length,gene_count,PHYLUM,SPECIES)) %>% arrange(coding_density)
write_csv(highdensity,file.path(statsplotdir,"highdensity_chroms.csv"))
# closeup shop

CBM18_sql = "SELECT sp.SPECIES, sp.LOCUSTAG, pfam_id, COUNT(*) as domain_count
FROM species as sp, pfam
WHERE pfam.species_prefix = sp.LOCUSTAG AND
pfam.pfam_id = 'Chitin_bind_1'
GROUP BY sp.SPECIES, sp.LOCUSTAG, pfam_id ORDER by domain_count DESC"

CBM18_res <- dbGetQuery(con, CBM18_sql)
head(CBM18_res)

dbDisconnect(con, shutdown = TRUE)


