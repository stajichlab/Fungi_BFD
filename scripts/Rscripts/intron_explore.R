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

# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir="intronDB/introns.duckdb", read_only = TRUE)

speciessql = "SELECT s.LOCUSTAG, PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN, TOTAL_LENGTH, N50
FROM species as s, asm_stats as stats
WHERE s.LOCUSTAG = stats.LOCUSTAG"

speciesinfo <- dbGetQuery(con, speciessql)

genelensql="
SELECT g.*, PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN
FROM species s,
(SELECT LOCUSTAG, avg(length) as mean_protein_length, median(length) as median_protein_length, 
                 avg(cds_length) as mean_cds_length, median(cds_length) as median_cds_length
FROM (SELECT substring(transcript_id,1,8) as LOCUSTAG, p.length, p.length * 3 as cds_length FROM gene_proteins p)
GROUP BY LOCUSTAG) as g
WHERE g.LOCUSTAG = s.LOCUSTAG"
genelen <- dbGetQuery(con, genelensql)
head(genelen)
# (SELECT abs(gene_introns.end - gene_introns.start) as intron_length

# rewrite to add introns per KB 
introncountsql="
SELECT i.*, PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN
FROM 
(SELECT LOCUSTAG, avg(intron_count) as mean_intron_ct, 
        median(intron_count) as median_intron_ct
 FROM (SELECT substring(transcript_id,1,8) as LOCUSTAG, 
              max(intron_number + 1) as intron_count, transcript_id
      FROM gene_introns 
      GROUP BY transcript_id)
 GROUP BY LOCUSTAG) as i,
species s
WHERE
i.LOCUSTAG = s.LOCUSTAG
"

intronct <- dbGetQuery(con, introncountsql)
head(intronct)
## test
testsql = "SELECT max(gene_introns.intron_number + 1) as intron_count, gene_introns.transcript_id
      FROM gene_introns WHERE transcript_id LIKE 'FE0E32D4_006065%'
      GROUP BY gene_introns.transcript_id"
testdat <- dbGetQuery(con, testsql)
testdat
testlensql = "SELECT transcript_id, p.length, p.length * 3 as cds_length FROM gene_proteins p
              WHERE transcript_id LIKE 'FE0E32D4_006065%'"
testlendat <- dbGetQuery(con, testlensql)
testlendat
### 

# intron length calculated
intronlensql="
SELECT gi.*, PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN
FROM species s,
 (SELECT LOCUSTAG, avg(intron_length) as mean_intron_len, median(intron_length) as median_intron_len
  FROM (SELECT substring(transcript_id,1,8) as LOCUSTAG, 
               abs(gene_introns.end - gene_introns.start) as intron_length
        FROM gene_introns)
  GROUP BY LOCUSTAG) as gi
WHERE gi.LOCUSTAG = s.LOCUSTAG"

intronlens <- dbGetQuery(con, intronlensql)
head(intronlens)

# intron frequency calculated
intronfreqsql="
SELECT gpk.*, PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN
FROM species s,
 ( SELECT LOCUSTAG, avg(intronsperkb) as mean_intronsperkb,
          median(intronsperkb) as median_intronsperkb,
   FROM (SELECT introns.LOCUSTAG, p.transcript_id, 1000 * IFNULL(intron_count,0) / (3 * p.length) as intronsperkb
         FROM gene_proteins AS p LEFT JOIN 
              (SELECT substring(gene_introns.transcript_id,1,8) as LOCUSTAG, 
                      max(gene_introns.intron_number + 1) as intron_count,
                      gene_introns.transcript_id
              FROM gene_introns 
              GROUP BY gene_introns.transcript_id) as introns 
         ON p.transcript_id = introns.transcript_id)
  GROUP BY LOCUSTAG) as gpk
WHERE gpk.LOCUSTAG = s.LOCUSTAG"

intronfreq <- dbGetQuery(con, intronfreqsql)
head(intronfreq)

intronfreq$SUBPHYLUM = factor(intronfreq$SUBPHYLUM)
intronfreq$PHYLUM = factor(intronfreq$PHYLUM)
intronfreq$SUBPHYLUM = factor(intronfreq$SUBPHYLUM,
        levels = as.character(intronfreq$SUBPHYLUM)[order(intronfreq$PHYLUM)])

sumbysubphylum_freq <- intronfreq %>%
  filter(SUBPHYLUM != "Ascomycotina" & ! is.na(SUBPHYLUM)) %>%  # remove this is too small numbers
  group_by(SUBPHYLUM,PHYLUM) %>%
  summarise(freq_mean = mean(mean_intronsperkb),
            freq_sd = sd(mean_intronsperkb),
            freq_N = n(),
            freq_se = freq_sd / sqrt(freq_N),
            freq_upper_limit = freq_mean + freq_se,
            freq_lower_limit = freq_mean - freq_se,
            freq_total = sum(mean_intronsperkb),
            freq_median = median(median_intronsperkb)) %>%
  filter( freq_N >= 3)


intronlensql="
SELECT il.*, PHYLUM, SUBPHYLUM, CLASS, s.ORDER, GENUS, s.SPECIES, s.STRAIN, 
      asm_stats.TOTAL_LENGTH, asm_stats.CONTIG_COUNT
FROM species s, asm_stats,
 (SELECT LOCUSTAG, avg(intron_length) as mean_length, median(intron_length) as median_length
 FROM (SELECT substring(transcript_id,1,8) as LOCUSTAG, 
        (gene_introns.end - gene_introns.start) as intron_length
      FROM gene_introns) 
  GROUP BY LOCUSTAG) as il
WHERE s.LOCUSTAG = il.LOCUSTAG and asm_stats.LOCUSTAG = s.LOCUSTAG
"

intronlen <- dbGetQuery(con, intronlensql)
head(intronlen)

intronlen_sum <- intronlen %>% 
  filter(SUBPHYLUM != "Ascomycotina" & ! is.na(SUBPHYLUM)) %>%  # remove this is too small numbers  
  group_by(SUBPHYLUM, PHYLUM) %>% 
   summarise(ilen_mean = mean(mean_length),
             ilen_sd = sd(mean_length),
             ilen_N = n(),
             ilen_se = ilen_sd / sqrt(ilen_N),
             ilen_upper_limit = ilen_mean + ilen_se,
             ilen_lower_limit = ilen_mean - ilen_se,
             ilen_total = sum(mean_length),
             ilen_median = median(mean_length)) %>% filter(ilen_N >= 3)
   


subphylumdata <- sumbysubphylum_freq %>% inner_join(intronlen_sum %>% 
                                                select(-c(PHYLUM)),by="SUBPHYLUM")

subphylumdata$SUBPHYLUM = factor(subphylumdata$SUBPHYLUM)
subphylumdata$PHYLUM = factor(subphylumdata$PHYLUM, 
                              ordered=TRUE,
                              levels=unique(sort(subphylumdata$PHYLUM)))

subphylumdata$SUBPHYLUM = factor(subphylumdata$SUBPHYLUM,
                                 ordered=TRUE,
                                 levels = as.character(subphylumdata$SUBPHYLUM)[order(subphylumdata$PHYLUM)])

subphylumdata <- subphylumdata %>%
  group_by(PHYLUM) %>%
  mutate(lvl = as.numeric(factor(SUBPHYLUM, levels = unique(SUBPHYLUM)))) %>%
  ungroup()


# Define the number of colors you want
nb.cols <- 25
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

plotdata <- subphylumdata %>% mutate(main_color = brewer_pal("qual")(n_distinct(subphylumdata$PHYLUM))[PHYLUM] ) %>%
  arrange(PHYLUM, lvl) %>%
  mutate( relvl = rescale(lvl, c(0, .6), 
                          from = c(0, max(length(lvl), 4))),
          sub_color= darken(main_color, relvl))


p <- ggplot(plotdata,aes(x=freq_mean,
                         y=ilen_mean,
                         color=sub_color,
                         fill=sub_color)) +
  geom_point() +
  geom_label_repel(aes(label = SUBPHYLUM),
                   color="black",
                   size = 4,
                   nudge_x = -0.2, direction = "y", hjust = "right",
                   box.padding   = 0.1, 
                   point.padding = 0.1,
                   max.overlaps = 15,
                   segment.color = 'grey50') +
  geom_errorbar(
    aes(ymin=ilen_lower_limit, 
        ymax=ilen_upper_limit,
        xmin=freq_lower_limit,
        xmax=freq_upper_limit,
      ), 
    width=.2) +
  scale_fill_identity() +
  scale_color_identity() +
  xlab("Mean introns per kb") +
  ylab("Mean intron length") +
  theme_cowplot(12) 
p

ggsave("plots/intron_size_freq.pdf",p,width=15,height=10)

intronlen$PHYLUM <- factor(intronlen$PHYLUM)
intronlen$SUBPHYLUM <- factor(intronlen$SUBPHYLUM)

p2 <- ggplot(intronlen,aes(x=TOTAL_LENGTH,
                         y=mean_length,
                         color=SUBPHYLUM,
                         shape=PHYLUM)) +
  scale_shape_manual(values=1:nlevels(intronlen$PHYLUM)) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  geom_point() +
  xlab("Genome Size") +
  ylab("Mean intron length") + scale_x_log10() +
  theme_cowplot(12) 
p2

ggsave("plots/intronlen_vs_genomesize.pdf",p2,width=10,height=10)

sumbysubphylum_count <- intronct %>% left_join(speciesinfo,by="LOCUSTAG") %>%
  group_by(SUBPHYLUM) %>% summarise(ct_mean = mean(intron_count),
                                    ct_sd = sd(intron_count),
                                    ct_n = n(),
                                    ct_total = sum(intron_count),
                                    ct_median = median(intron_count))


dbDisconnect(con, shutdown = TRUE)


