#!/usr/bin/env Rscript

library(tidyverse)
library(duckdb)
library(purrr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)

plotdir = "plots/function_explore"
dir.create(file.path(plotdir),showWarnings = FALSE)
pdf(file.path(plotdir,"IDR_domain_misc.pdf"))

# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir="functionalDB/function.duckdb", read_only = TRUE)

asmstats_sql = "SELECT sp.SUBPHYLUM, sp.GENUS, asm_stats.*
FROM asm_stats, species sp
WHERE
sp.LOCUSTAG = asm_stats.LOCUSTAG"
asm_stats <- dbGetQuery(con, asmstats_sql)

IDR_sql = "select genecount.LOCUSTAG, species.PHYLUM, species.SUBPHYLUM, species.GENUS, species.SPECIES, asm_stats.TOTAL_LENGTH, asm_stats.GC_PERCENT, IDP.idp_count, 
           genecount.pep_count, 100 * IDP.idp_count/genecount.pep_count as percent_IDP 
FROM species, asm_stats,
     (select LOCUSTAG, count(*) as pep_count FROM gene_proteins GROUP BY LOCUSTAG) as genecount, 
     (select species_prefix as sp, count(*) as idp_count from idp_summary where IDP_fraction > 0.8 GROUP BY species_prefix) as IDP
WHERE genecount.LOCUSTAG=IDP.sp AND genecount.LOCUSTAG = species.LOCUSTAG and asm_stats.LOCUSTAG = species.LOCUSTAG
ORDER BY percent_IDP DESC"
IDP_sumstat <- dbGetQuery(con, IDR_sql)
IDP_sumstat
nb.cols <- length(unique(IDP_sumstat$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

p <- ggplot(IDP_sumstat) +
  geom_point(aes(x=TOTAL_LENGTH/1000000, y=idp_count, color=SUBPHYLUM, fill=SUBPHYLUM),size=2) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("IDP count") +
  xlab("Genome size (Mb)") + ggtitle("IDP vs Genome Size") +
  theme_cowplot(12) + theme(legend.position="bottom")
p

p <- ggplot(IDP_sumstat) +
  geom_point(aes(x=pep_count, y=idp_count, color=SUBPHYLUM, fill=SUBPHYLUM),size=2) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("IDP count") +
  xlab("Protein count") + ggtitle("IDP vs Protein count") +
  theme_cowplot(12) + theme(legend.position="bottom")
p

sumval <- IDP_sumstat %>% group_by(SUBPHYLUM) %>% summarise(
  sp_mean = mean(percent_IDP),
  sp_N = n(),
  sp_sd = sd(percent_IDP),
  sp_se = sp_sd / sqrt(sp_N))

p <- ggplot(IDP_sumstat) +
  geom_bar(data=sumval,
           aes(y=sp_mean,
               x=SUBPHYLUM,
               fill=SUBPHYLUM,
               color=SUBPHYLUM),
           stat="identity", width=0.5) +
  geom_point(aes(x=SUBPHYLUM, y=percent_IDP, color=SUBPHYLUM),size=2) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("IDP Frequency") +
  xlab("Subphylum") + ggtitle("IDP frequency per Subphylum") +
  theme_cowplot(12) + theme(legend.position="bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
mean(IDP_sumstat$percent_IDP)

ggsave(file.path(plotdir,"IDP_summary.pdf"),p,width=10,height=8)
dbDisconnect(con, shutdown = TRUE)


dbDisconnect(con, shutdown = TRUE)


