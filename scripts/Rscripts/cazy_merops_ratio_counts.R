#!/usr/bin/env Rscript

library(tidyverse)
library(duckdb)
library(purrr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
pdf("plots/cazy_plots_misc.pdf")
# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir="functionalDB/function.duckdb", read_only = TRUE)

asmstats_sql = "SELECT sp.SUBPHYLUM, sp.GENUS, asm_stats.*
FROM asm_stats, species sp
WHERE
sp.LOCUSTAG = asm_stats.LOCUSTAG"
asm_stats <- dbGetQuery(con, asmstats_sql)

# barplot for these GENERA

sumsizes <- asm_stats %>% group_by(GENUS) %>% summarise(
  g_mean = mean(TOTAL_LENGTH / 1000000),
  g_N = n(),
  g_sd = sd(TOTAL_LENGTH/1000000),
  g_se = g_sd / sqrt(g_N))

speciessql = "SELECT s.LOCUSTAG, PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN, TOTAL_LENGTH, N50
FROM species as s, asm_stats as stats
WHERE s.LOCUSTAG = stats.LOCUSTAG"

speciesinfo <- dbGetQuery(con, speciessql)

# PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN,
# species as s, 
# WHERE m.species_prefix = s.LOCUSTAG
cazymerops_sql="
SELECT s.LOCUSTAG, m.merops_count, c.cazy_count, m.merops_count/c.cazy_count as MERCAZ_ratio, PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN,
FROM species s,
(SELECT species_prefix, COUNT(*) as merops_count,
FROM (SELECT DISTINCT species_prefix, protein_id FROM merops) as m
GROUP BY m.species_prefix) as m,
(SELECT species_prefix, COUNT(*) as cazy_count,
FROM (SELECT DISTINCT species_prefix, protein_id FROM cazy) as c
GROUP BY c.species_prefix) as c
WHERE m.species_prefix = s.LOCUSTAG and c.species_prefix = s.LOCUSTAG"
cazymerops <- dbGetQuery(con, cazymerops_sql)

head(cazymerops)

CM_subphylum <- cazymerops %>% group_by(SUBPHYLUM,PHYLUM) %>% 
  summarize(cm_mean = mean(MERCAZ_ratio),
            cm_sd = sd(MERCAZ_ratio),
            cm_N = n(),
            cm_se = cm_sd / sqrt(cm_N),
            cm_upper_limit = cm_mean + cm_se,
            cm_lower_limit = cm_mean - cm_se,
            cm_total = sum(MERCAZ_ratio),
            cm_median = median(MERCAZ_ratio))


CM_subphylum$SUBPHYLUM = factor(CM_subphylum$SUBPHYLUM)
CM_subphylum$PHYLUM = factor(CM_subphylum$PHYLUM)
CM_subphylum$SUBPHYLUM = factor(CM_subphylum$SUBPHYLUM,
                                        levels = as.character(CM_subphylum$SUBPHYLUM)[order(CM_subphylum$PHYLUM)])
# Define the number of colors you want
nb.cols <- 9
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

p <- ggplot() + geom_bar(data=CM_subphylum,
                         aes(y=cm_mean,
                             x=SUBPHYLUM,
                             fill=PHYLUM,
                             color=PHYLUM),
                         stat="identity", width=0.5) +
  geom_errorbar(data = CM_subphylum,
                aes(y=cm_mean,
                    x=SUBPHYLUM,
                    ymin = cm_mean - cm_se,
                    ymax = cm_mean + cm_se,
                    color = PHYLUM),
                stat="identity", width=0.5) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  xlab("Subphylum") +
  ylab("MEROPS/CAZY ratio [all]") +
  theme_cowplot(12) + 
  ggtitle("MEROPS / CAZY ratio by Subphylum") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14))


#p
ggsave("plots/MEROPS_CAZY_profile.pdf",p,width=10,height=8)

cazymeropssecreted_sql="
SELECT s.LOCUSTAG, m.merops_count, c.cazy_count, m.merops_count/c.cazy_count as MERCAZ_ratio,
PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN,
FROM species s,
(SELECT m.species_prefix, COUNT(*) as merops_count
FROM (SELECT DISTINCT species_prefix, protein_id FROM merops) as m, signalp sp
WHERE sp.protein_id = m.protein_id
GROUP BY (m.species_prefix,sp.species_prefix)) as m,
(SELECT sp.species_prefix, COUNT(*) as cazy_count
FROM (SELECT DISTINCT species_prefix, protein_id FROM cazy) as c, signalp sp
WHERE sp.protein_id = c.protein_id
GROUP BY (c.species_prefix,sp.species_prefix)) as c
WHERE m.species_prefix = s.LOCUSTAG and c.species_prefix = s.LOCUSTAG"
cazymeropssecreted <- dbGetQuery(con, cazymeropssecreted_sql)
head(cazymeropssecreted)

CMsecreted_subphylum <- cazymeropssecreted %>% group_by(SUBPHYLUM,PHYLUM) %>% 
  summarize(cm_mean = mean(MERCAZ_ratio),
            cm_sd = sd(MERCAZ_ratio),
            cm_N = n(),
            cm_se = cm_sd / sqrt(cm_N),
            cm_upper_limit = cm_mean + cm_se,
            cm_lower_limit = cm_mean - cm_se,
            cm_total = sum(MERCAZ_ratio),
            cm_median = median(MERCAZ_ratio))


CMsecreted_subphylum$SUBPHYLUM = factor(CMsecreted_subphylum$SUBPHYLUM)
CMsecreted_subphylum$PHYLUM = factor(CMsecreted_subphylum$PHYLUM)
CMsecreted_subphylum$SUBPHYLUM = factor(CMsecreted_subphylum$SUBPHYLUM,
                                        levels = as.character(CMsecreted_subphylum$SUBPHYLUM)[order(CMsecreted_subphylum$PHYLUM)])
# Define the number of colors you want
nb.cols <- 9
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)


p <- ggplot() + geom_bar(data=CMsecreted_subphylum,
                         aes(y=cm_mean,
                             x=SUBPHYLUM,
                             fill=PHYLUM,
                             color=PHYLUM),
                         stat="identity", width=0.5) +
  geom_errorbar(data = CMsecreted_subphylum,
                aes(y=cm_mean,
                    x=SUBPHYLUM,
                    ymin = cm_mean - cm_se,
                    ymax = cm_mean + cm_se,
                    color = PHYLUM),
                stat="identity", width=0.5) +
  #geom_point(data=cazymeropssecreted, aes(y=MERCAZ_ratio, x=SUBPHYLUM, fill=SUBPHYLUM),size=2) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  xlab("Subphylum") +
  ylab("MEROPS/CAZY [secreted only] ratio") +
  theme_cowplot(12) + 
  ggtitle("MEROPS / CAZY ratio for proteins with SignalP signal") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14))
  

p
ggsave("plots/MEROPS_CAZY_profile_secreted.pdf",p,width=10,height=8)

# geom_point(data=asm_stats, aes(y=TOTAL_LENGTH/1000000, x=SUBPHYLUM, fill=SUBPHYLUM),size=2) +

nb.cols <- length(unique(cazymerops$PHYLUM))
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)


p <- ggplot(cazymerops) + geom_point(aes(x=cazy_count,
                             y=merops_count,
                             fill=PHYLUM,
                             color=PHYLUM), alpha=0.75) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("MEROPS Count") +
  xlab("CAZY Count") +
  theme_cowplot(12) + 
  ggtitle("MEROPS vs CAZY counys for all proteins") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14))


p

digestivequery_sql = 
"
SELECT ssp.LOCUSTAG, ssp.ssp_count, secreted.secreted_count, de.de_count, merops.merops_count, funguild.trophicMode, species.PHYLUM, species.SUBPHYLUM, species.SPECIES, genes.gene_count,
       asm_stats.TOTAL_LENGTH
FROM 

(SELECT gene_info.LOCUSTAG, COUNT(DISTINCT gene_info.gene_id) AS secreted_count
FROM signalp, gene_info, gene_proteins
WHERE signalp.protein_id = gene_info.gene_id AND gene_proteins.gene_id = gene_info.gene_id AND
signalp.probability > 0.60 AND signalp.protein_id not in (select cazy.protein_id FROM cazy where cazy.coverage > 0.50 AND cazy.evalue < 1e-5 AND cazy.HMM_id NOT LIKE 'GT%' AND cazy.HMM_id NOT LIKE 'fungi_doc%')
GROUP BY gene_info.LOCUSTAG) as secreted,

(SELECT gene_info.LOCUSTAG, COUNT(DISTINCT gene_info.gene_id) AS ssp_count
FROM signalp, gene_info, gene_proteins
WHERE signalp.protein_id = gene_info.gene_id AND gene_proteins.gene_id = gene_info.gene_id AND
signalp.probability > 0.60 AND gene_proteins.length < 300 
AND signalp.protein_id not in (select cazy.protein_id FROM cazy where cazy.coverage > 0.50 AND cazy.evalue < 1e-5 AND cazy.HMM_id NOT LIKE 'GT%' AND cazy.HMM_id NOT LIKE 'fungi_doc%')
GROUP BY gene_info.LOCUSTAG) as ssp,

(SELECT gene_info.LOCUSTAG, COUNT(DISTINCT gene_info.gene_id) AS de_count
FROM cazy, gene_info, gene_proteins
WHERE cazy.protein_id = gene_info.gene_id AND gene_proteins.gene_id = gene_info.gene_id AND
cazy.coverage > 0.50 AND cazy.evalue < 1e-5 AND cazy.HMM_id NOT LIKE 'GT%' AND cazy.HMM_id NOT LIKE 'fungi_doc%'
GROUP BY gene_info.LOCUSTAG) as de,

(SELECT gene_info.LOCUSTAG, COUNT(DISTINCT gene_info.gene_id) AS merops_count
FROM merops, gene_info, gene_proteins
WHERE merops.protein_id = gene_info.gene_id AND gene_proteins.gene_id = gene_info.gene_id AND
merops.aln_length / gene_proteins.length > 0.50 AND merops.evalue < 1e-10
GROUP BY gene_info.LOCUSTAG) as merops,

(SELECT gene_info.LOCUSTAG, COUNT(DISTINCT gene_info.gene_id) AS gene_count
FROM gene_info
GROUP BY gene_info.LOCUSTAG) as genes,

funguild, species, asm_stats

WHERE ssp.LOCUSTAG = de.LOCUSTAG AND ssp.LOCUSTAG = merops.LOCUSTAG AND 
ssp.LOCUSTAG = funguild.species_prefix AND species.LOCUSTAG = ssp.LOCUSTAG AND
genes.LOCUSTAG = ssp.LOCUSTAG AND asm_stats.LOCUSTAG = genes.LOCUSTAG and secreted.LOCUSTAG = ssp.LOCUSTAG
"



digestive <- dbGetQuery(con, digestivequery_sql)
head(digestive)
digestive %>% filter(LOCUSTAG == "Podan3")

zygoonly = digestive %>% filter(PHYLUM == "Zoopagomycota" | PHYLUM == "Mucoromycota")
nb.cols <- length(unique(zygoonly$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)


p <- ggplot(zygoonly) + geom_point(aes(x=log(de_count / (secreted_count - de_count))/log(10),
                                         y=log(ssp_count / (secreted_count - ssp_count))/log(10),
                                         fill=SUBPHYLUM,
                                         color=SUBPHYLUM,
                                         size=gene_count,
                                         ), alpha=0.75) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("SSP Count/(Total - SSP Count)") +
  xlab("DE Count/(Total - DE Count)") +
  theme_cowplot(12) + 
  ggtitle("SSP vs CAZY counts, Mucoromycota and Zoopagomycota") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14))


p
ggsave("plots/1kfg_zygo_SSP_vs_DE_gene_count.pdf",p,width=10,height=8)
nb.cols <- length(unique(digestive$PHYLUM))
mycolors <- colorRampPalette(brewer.pal(9, "Set1"))(nb.cols)

p <- ggplot(digestive) + geom_point(aes(x=log(de_count / (secreted_count - de_count))/log(10),
                                        y=log(ssp_count / (secreted_count - ssp_count))/log(10),
                                        fill=PHYLUM,
                                        color=PHYLUM,
                                        size=gene_count,
), alpha=0.75) +
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("SSP Count/(Total - SSP Count)") +
  xlab("DE Count/(Total - DE Count)") +
  theme_cowplot(12) + 
  ggtitle("SSP vs CAZY counts for all proteins") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14))


p
ggsave("plots/1kfg_all_SSP_vs_DE_gene_count.pdf",p,width=10,height=8)
dbDisconnect(con, shutdown = TRUE)


dbDisconnect(con, shutdown = TRUE)


