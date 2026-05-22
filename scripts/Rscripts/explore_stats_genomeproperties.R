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
statsplotdir = file.path(plotdir,"asm_seqstats")
dir.create(statsplotdir, showWarnings = FALSE)

pdf(file.path(statsplotdir,"plot_misc.pdf"))
# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir=file.path(DBDIR,DBNAME), read_only = TRUE)

asmstat_sql ="
SELECT sp.*, stats.GC_PERCENT, stats.TOTAL_LENGTH, gene_count, mean_gene_length
FROM 
species as sp,
asm_stats as stats,
(SELECT LOCUSTAG, count(*) as gene_count, MEAN(length) as mean_gene_length
FROM gene_proteins 
GROUP BY LOCUSTAG) as gp
WHERE 
sp.LOCUSTAG = stats.LOCUSTAG and gp.LOCUSTAG = sp.LOCUSTAG"


asmstat_res <- dbGetQuery(con, asmstat_sql)
head(asmstat_res)

# barplot for these GENERA
sumAsm = asmstat_res %>% group_by(ORDER) %>% 
  summarize(
  gc_mean = mean(GC_PERCENT),
  gc_N = n(),
  gc_sd = sd(GC_PERCENT),
  gc_se = gc_sd / sqrt(gc_N),
  len_mean = mean(TOTAL_LENGTH / 1000000),
  len_N = n(),
  len_sd = sd(TOTAL_LENGTH/1000000),
  len_se = len_sd / sqrt(len_N))

len_p <- ggplot(sumAsm) + geom_bar(
                         aes(y=len_mean,
                             x=ORDER),
                         stat="identity", width=0.5) +
  geom_errorbar(
                aes(y=len_mean,
                    x=ORDER,
                    ymin = len_mean - len_se,
                    ymax = len_mean + len_se),
                stat="identity", width=0.5) +
  geom_point(aes(y=len_mean, x=ORDER),size=2) + 
  xlab("ORDER") +
  ylab("Genome size (Mb)") +
  theme_cowplot(12) + theme(axis.text.x = element_text(angle = 90,size = 6))

len_p

ggsave(file.path(statsplotdir,"genome_size_order.pdf"),len_p,width=24,height=10)

nb.cols <- length(unique(asmstat_res$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

count_glen_p <- ggplot(asmstat_res) +
  geom_point(aes(x=gene_count, y=mean_gene_length,color=SUBPHYLUM, fill=SUBPHYLUM),size=1.5,alpha=0.5) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle("Gene Count vs Mean Gene length") +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10() 
count_glen_p
ggsave(file.path(statsplotdir,"genecount_gene_length.pdf"),count_glen_p,width=15,height=10)

phylum_countglen_p <- ggplot(asmstat_res) +
  geom_point(aes(x=gene_count, y=mean_gene_length,color=SUBPHYLUM, fill=SUBPHYLUM),size=1.5,alpha=0.5) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle("Gene Count vs Mean Gene length") +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10() + facet_wrap(~PHYLUM)

phylum_countglen_p
ggsave(file.path(statsplotdir,"genecount_gene_length_facet.pdf"),phylum_countglen_p,width=15,height=10)

basidio <- asmstat_res %>% filter(PHYLUM=="Basidiomycota")
nb.cols <- length(unique(basidio$CLASS))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

basidio_countglen_p <- ggplot(basidio) +
  geom_point(aes(x=gene_count, y=mean_gene_length,color=CLASS, fill=CLASS),size=1.5,alpha=0.7) + 
  geom_smooth(method = "lm", se = FALSE,color="black",aes(x=gene_count, y=mean_gene_length),formula = y ~ x) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle("Gene Count vs Mean Gene length") +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10() + facet_wrap(~SUBPHYLUM)

basidio_countglen_p
ggsave(file.path(statsplotdir,"genecount_gene_length_basidio.pdf"),basidio_countglen_p,width=15,height=10)

# could make this a function + lapply
for (subphylum in unique(basidio$SUBPHYLUM))
{
  subph <- basidio %>% filter(SUBPHYLUM==subphylum)

  countglen_p <- ggplot(subph,aes(x=gene_count, y=mean_gene_length)) +
  geom_point(aes(color=CLASS, fill=CLASS),size=1.5,alpha=0.7) + 
    geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
    scale_colour_brewer(palette = "Dark2") +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle(sprintf("%s: Gene Count vs Mean Gene length",subphylum)) +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10()

  countglen_p
  ggsave(file.path(statsplotdir,sprintf("genecount_gene_length_basidio_%s.pdf",subphylum)),countglen_p,width=8,height=8)
}

mucoro <- asmstat_res %>% filter(PHYLUM=="Mucoromycota")

  countglen_p <- ggplot(mucoro,aes(x=gene_count, y=mean_gene_length)) +
    geom_point(aes(color=CLASS, fill=CLASS,shape=SUBPHYLUM),size=2,alpha=0.7) + 
    geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
    scale_colour_brewer(palette = "Set2") +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle(sprintf("%s: Gene Count vs Mean Gene length",unique(mucoro$PHYLUM))) +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10()

  countglen_p
  ggsave(file.path(statsplotdir,"genecount_gene_length_mucoromycota.pdf"),countglen_p,width=12,height=8)

zoopag <- asmstat_res %>% filter(PHYLUM=="Zoopagomycota")
# switch this around and just plot all Zooags together
# could make this a function + lapply
  countglen_p <- ggplot(zoopag %>% filter(!is.na(CLASS)),aes(x=gene_count, y=mean_gene_length)) +
  geom_point(aes(color=CLASS, fill=CLASS,shape=SUBPHYLUM),size=2,alpha=0.7) + 
    geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
    scale_colour_brewer(palette = "Set2") +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle(sprintf("%s: Gene Count vs Mean Gene length",unique(zoopag$PHYLUM))) +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10()

  countglen_p
  ggsave(file.path(statsplotdir,"genecount_gene_length_zoopag.pdf"),countglen_p,width=8,height=8)


 chytrid <- asmstat_res %>% filter(PHYLUM=="Chytridiomycota")
  countglen_p <- ggplot(chytrid,aes(x=gene_count, y=mean_gene_length)) +
    geom_point(aes(color=CLASS, fill=CLASS,shape=SUBPHYLUM),size=2,alpha=0.7) + 
    geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
    scale_colour_brewer(palette = "Set2") +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle(sprintf("%s: Gene Count vs Mean Gene length",unique(chytrid$PHYLUM))) +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10()

  countglen_p
  ggsave(file.path(statsplotdir,"genecount_gene_length_chytrid.pdf"),countglen_p,width=8,height=8)

asco <- asmstat_res %>% filter(PHYLUM=="Ascomycota")
for (subphylum in unique(asco$SUBPHYLUM))
{
  subph <- asco %>% filter(SUBPHYLUM==subphylum)
  nb1.cols <- length(unique(subph$CLASS))
  mycolors1 <- colorRampPalette(brewer.pal(8, "Set1"))(nb1.cols)
  countglen_p <- ggplot(subph,aes(x=gene_count, y=mean_gene_length)) +
    geom_point(aes(color=CLASS, fill=CLASS),size=1.5,alpha=0.7) + 
    geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
    scale_colour_manual(values = mycolors1) +
    scale_fill_manual(values = mycolors1) +
    ylab("Gene Length") +
    xlab("Gene count (Mb)") + ggtitle(sprintf("%s: Gene Count vs Mean Gene length",subphylum)) +
    theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10()
  
  countglen_p
  ggsave(file.path(statsplotdir,sprintf("genecount_gene_length_asco_%s.pdf",subphylum)),countglen_p,width=8,height=8)
}

asco <- asmstat_res %>% filter(PHYLUM=="Ascomycota") %>% filter(SUBPHYLUM != "NA")
countglen_p <- ggplot(asco,aes(x=gene_count, y=mean_gene_length)) +
  geom_point(aes(color=SUBPHYLUM, fill=SUBPHYLUM),size=1.5,alpha=0.7) + 
  geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
  scale_colour_brewer(palette = "Set1") + 
    ylab("Gene Length") +
    xlab("Gene count (Mb)") + ggtitle(sprintf("%s: Gene Count vs Mean Gene length","Ascomycota")) +
    theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10()
  
  countglen_p
ggsave(file.path(statsplotdir,"genecount_gene_length_asco_all.pdf"),countglen_p,width=8,height=8)

# could make this a function + lapply
for (subphylum in unique(basidio$SUBPHYLUM))
{
  subph <- basidio %>% filter(SUBPHYLUM==subphylum)
  countglen_p <- ggplot(subph,aes(x=gene_count, y=mean_gene_length)) +
  geom_point(aes(color=CLASS, fill=CLASS),size=1.5,alpha=0.7) + 
    geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
    scale_colour_brewer(palette = "Dark2") +
  ylab("Gene Length") +
  xlab("Gene count (Mb)") + ggtitle(sprintf("%s: Gene Count vs Mean Gene length",subphylum)) +
  theme_cowplot(12) + theme(legend.position="bottom") + scale_x_log10()

  countglen_p
  ggsave(file.path(statsplotdir,sprintf("genecount_gene_length_basidio_%s.pdf",subphylum)),countglen_p,width=8,height=8)
}

## Codon Usage plots

codonfreq_sql ="
SELECT sp.*, cf.*, stats.GC_PERCENT, stats.TOTAL_LENGTH
FROM 
species as sp,
codon_frequency as cf,
asm_stats as stats
WHERE 
sp.LOCUSTAG = cf.species_prefix AND
sp.LOCUSTAG = stats.LOCUSTAG"

aafreq_sql ="
SELECT sp.*, aaf.*, stats.GC_PERCENT, stats.TOTAL_LENGTH
FROM 
species as sp,
aa_frequency as aaf,
asm_stats as stats
WHERE 
sp.LOCUSTAG = aaf.species_prefix AND
sp.LOCUSTAG = stats.LOCUSTAG"

codonfreq_res <- dbGetQuery(con, codonfreq_sql)
head(codonfreq_res)

aafreq_res <- dbGetQuery(con, aafreq_sql)
head(aafreq_res)
aafreq_wide <- aafreq_res %>% select(c(PHYLUM,SUBPHYLUM,CLASS,GENUS,SPECIES,GC_PERCENT,TOTAL_LENGTH,LOCUSTAG,amino_acid, frequency)) %>% 
  filter(! is.na(PHYLUM) ) %>% filter( ! PHYLUM %in% c("Sanchytriomycota","Cryptomycota")) %>%
  pivot_wider(id_cols = c(PHYLUM,SUBPHYLUM,CLASS,GENUS,SPECIES,GC_PERCENT,TOTAL_LENGTH,LOCUSTAG), 
              names_from = amino_acid, values_from = frequency)
aa_pcadat <- as.matrix(aafreq_wide %>% select(-c(PHYLUM,SUBPHYLUM,CLASS,GENUS,SPECIES,GC_PERCENT,TOTAL_LENGTH,LOCUSTAG,X)))
rownames(aa_pcadat) <- aafreq_wide$LOCUSTAG
aa_pca_res <- prcomp(aa_pcadat, scale. = TRUE)

nb.cols <- length(unique(aafreq_wide$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

aa_pcaplot<- autoplot(aa_pca_res, data = aafreq_wide, shape='PHYLUM',colour = 'SUBPHYLUM', alpha=0.7,
                      label = FALSE, label.size = 3) + 
  scale_shape_manual(values=seq(0,8)) +
  theme_cowplot(12) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors)
#  scale_colour_brewer(palette = "Set1") +
aa_pcaplot
ggsave(file.path(statsplotdir,"PCA_aa_freq_all.pdf"),aa_pcaplot,width=14,height=14)

ascoonly <- aafreq_wide %>% filter(PHYLUM == "Ascomycota" & SUBPHYLUM != "NA")
aa_asco_pcadat <- as.matrix(ascoonly %>% select(-c(PHYLUM,SUBPHYLUM,CLASS,GENUS,SPECIES,GC_PERCENT,TOTAL_LENGTH,LOCUSTAG,X)))
rownames(aa_asco_pcadat) <- ascoonly$LOCUSTAG
aa_asco_pca_res <- prcomp(aa_asco_pcadat, scale. = TRUE)

nb.cols <- length(unique(ascoonly$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

aa_asco_pcaplot<- autoplot(aa_asco_pca_res, data = ascoonly, colour = 'SUBPHYLUM', alpha=0.7,
                      label = FALSE, label.size = 3) + 
  theme_cowplot(12) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) + 
  ggtitle("PCA plot of AA usage across Ascomycota")
aa_asco_pcaplot
ggsave(file.path(statsplotdir,"PCA_aa_freq_ASCOs.pdf"),aa_asco_pcaplot,width=14,height=14)

#### 
ascoonly <- aafreq_wide %>% filter(PHYLUM == "Ascomycota" & SUBPHYLUM != "NA")
aa_asco_pcadat <- as.matrix(ascoonly %>% select(-c(PHYLUM,SUBPHYLUM,CLASS,GENUS,SPECIES,GC_PERCENT,TOTAL_LENGTH,LOCUSTAG,X)))
rownames(aa_asco_pcadat) <- ascoonly$LOCUSTAG
aa_asco_pca_res <- prcomp(aa_asco_pcadat, scale. = TRUE)

#nb.cols <- length(unique(ascoonly$SUBPHYLUM))
#mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

aa_asco_pcaplot<- autoplot(aa_asco_pca_res, data = ascoonly, colour = 'SUBPHYLUM', alpha=0.7,
                      label = FALSE, label.size = 3) + 
  theme_cowplot(12) + 
  scale_colour_brewer(palette = "Set1") + 
  ggtitle("PCA plot of AA usage across Ascomycota")
aa_asco_pcaplot
ggsave(file.path(statsplotdir,"PCA_aa_freq_ASCOs.pdf"),aa_asco_pcaplot,width=14,height=14)

####
aa_pcafactors <- as_tibble(rownames_to_column(data.frame(aa_pca_res$x),var="LOCUSTAG"))

aa_pca_factors <- aafreq_wide %>% left_join(aa_pcafactors,by="LOCUSTAG")
fit <- lm(PC1~GC_PERCENT,aa_pca_factors%>%select(c(GC_PERCENT,PC1)))

aa_GC_plot <- ggplot(aa_pca_factors,aes(x=GC_PERCENT,y=PC1)) + geom_point(aes(color=PHYLUM,fill=PHYLUM)) + 
  theme_cowplot(12) + scale_colour_brewer(palette = "Set1") + 
  geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
  xlab("Genome GC %") +
  ylab("AA Freq PC1") + 
  ggtitle(paste("GC % vs AA Freq PC1",
                "Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                "Intercept =",signif(fit$coef[[1]],5 ),
                " Slope =",signif(fit$coef[[2]], 5),
                " p-value =",signif(summary(fit)$coef[2,4], 5))) +
  theme(legend.position="bottom") 
aa_GC_plot
ggsave(file.path(statsplotdir,"PCA_AA_freq_PC1_GC.pdf"),aa_GC_plot,width=14,height=14)

aa_pca_factors$SUBPHYLUM <- factor(aa_pca_factors$SUBPHYLUM)
aa_pca_factors$PHYLUM <- factor(aa_pca_factors$PHYLUM)
       

aa_GC_plot_f <- ggplot(aa_pca_factors,aes(x=GC_PERCENT,y=PC1)) + geom_point(aes(color=SUBPHYLUM)) + 
  theme_cowplot(12) + 
  geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
  xlab("Genome GC %") +
  ylab("AA Freq PC1") + 
  ggtitle("GC % vs AA Freq PC1") +
  theme(legend.position="bottom")  + facet_wrap(~PHYLUM )
aa_GC_plot_f
ggsave(file.path(statsplotdir,"PCA_AA_freq_PC1_GC_facet.pdf"),aa_GC_plot_f,width=14,height=14)
# CODON FREQ

codonfreq_wide <- codonfreq_res %>% select(c(PHYLUM,SUBPHYLUM,CLASS,GENUS,SPECIES,GC_PERCENT,TOTAL_LENGTH,LOCUSTAG,codon, frequency)) %>% 
  filter(! is.na(PHYLUM) ) %>% filter( ! PHYLUM %in% c("Sanchytriomycota","Cryptomycota")) %>%
  pivot_wider(id_cols = c(PHYLUM,SUBPHYLUM,CLASS,GENUS,SPECIES,GC_PERCENT,TOTAL_LENGTH,LOCUSTAG), 
              names_from = codon, values_from = frequency)
codon_pcadat <- as.matrix(codonfreq_wide %>% select(-c(PHYLUM,SUBPHYLUM,CLASS,GENUS,SPECIES,GC_PERCENT,TOTAL_LENGTH,LOCUSTAG)))
rownames(codon_pcadat) <- codonfreq_wide$LOCUSTAG
codon_pca_res <- prcomp(codon_pcadat, scale. = TRUE)

codon_pcaplot<- autoplot(codon_pca_res, data = codonfreq_wide, colour = 'PHYLUM', alpha=0.7, 
                         label = FALSE, label.size = 3) + 
  theme_cowplot(12) + scale_colour_brewer(palette = "Set1") 

ggsave(file.path(statsplotdir,"PCA_codon_freq_all.pdf"),codon_pcaplot,width=14,height=14)

codon_pcafactors <- as_tibble(rownames_to_column(data.frame(codon_pca_res$x),var="LOCUSTAG"))

codon_pca_factors <- codonfreq_wide %>% left_join(codon_pcafactors,by="LOCUSTAG")
fit <- lm(PC1~GC_PERCENT,codon_pca_factors%>%select(c(GC_PERCENT,PC1)))
cor2 <- cor(codon_pca_factors$GC_PERCENT,codon_pca_factors$PC1)

nb.cols <- length(unique(codon_pca_factors$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

codon_GC_plot <- ggplot(codon_pca_factors,aes(x=GC_PERCENT,y=PC1)) + geom_point(aes(color=SUBPHYLUM,fill=SUBPHYLUM)) + 
  theme_cowplot(12) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  #scale_colour_brewer(palette = "Set1") + 
  geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
  xlab("Genome GC %") +
  ylab("Codon Freq PC1") + 
  ggtitle(paste("GC % vs Codon Freq PC1",
    "Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                "Intercept =",signif(fit$coef[[1]],5 ),
                " Slope =",signif(fit$coef[[2]], 5),
                " p-value =",signif(summary(fit)$coef[2,4], 5))) +
  theme(legend.position="bottom") 
  
codon_GC_plot
ggsave(file.path(statsplotdir,"PCA_codon_freq_PC1_GC.pdf"),codon_GC_plot,width=14,height=14)


nb.cols <- length(unique(codon_pca_factors$SUBPHYLUM))
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)

codon_GC_plot_f <- ggplot(codon_pca_factors,aes(x=GC_PERCENT,y=PC1)) + 
  geom_point(aes(color=SUBPHYLUM)) + 
  scale_colour_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  theme_cowplot(12) + 
  geom_smooth(method = "lm", se = FALSE,color="black",formula = y ~ x) + 
  xlab("Genome GC %") +
  ylab("Codon Freq PC1") + 
  ggtitle("GC % vs Codon Freq PC1") +
  theme(legend.position="bottom")  + facet_wrap(~PHYLUM )
codon_GC_plot_f
ggsave(file.path(statsplotdir,"PCA_codon_freq_PC1_GC_facet.pdf"),codon_GC_plot_f,width=14,height=14)


# plot tRNA codon abundance vs codon usage
trnacodon_abundance_sql ="
SELECT sp.*, trna.codon as trna_codon, trna.trna_codon_count, trnatotal.trna_gene_count, 
       trna.trna_codon_count / trnatotal.trna_gene_count as relative_trna_abundance
FROM 
species sp,
(SELECT LOCUSTAG, codon, count(*) as trna_codon_count from gene_trna GROUP BY (LOCUSTAG,codon)) as trna,
(SELECT LOCUSTAG, count(*) as trna_gene_count from gene_trna GROUP BY (LOCUSTAG)) as trnatotal

WHERE 
sp.LOCUSTAG = trna.LOCUSTAG AND sp.LOCUSTAG = trnatotal.LOCUSTAG"

trnacodonfreq_res <- dbGetQuery(con, trnacodon_abundance_sql)
head(trnacodonfreq_res)

trna_codon_freq <- trnacodonfreq_res %>% select(c(LOCUSTAG,trna_codon,trna_codon_count,trna_gene_count,relative_trna_abundance)) %>% 
  left_join(codonfreq_res,by=c("trna_codon" = "codon",
                               "LOCUSTAG" = "species_prefix")) 
head(trna_codon_freq %>% arrange(LOCUSTAG,trna_codon))

# drop the no phylum data
trna_codon_freq <- trna_codon_freq %>% filter(!is.na(PHYLUM))

# frequency is codon frequency from CDS counts
# codoncount is the numner of 
trnacodonabun_p <- ggplot(trna_codon_freq %>% arrange(trna_codon, PHYLUM, SUBPHYLUM) %>%
                            filter(! is.na(frequency)),
                          aes(color=PHYLUM,fill=PHYLUM)) + 
  geom_point(aes(x=relative_trna_abundance,y=frequency)) + theme_cowplot(12)  + facet_wrap(~trna_codon,nrow=8,ncol=8) +
  scale_colour_brewer(palette = "Set1") + ylab("codon abundance") + xlab("relative tRNA abundance")

trnacodonabun_p

ggsave(file.path(statsplotdir,"tRNA_abundance_codon_frequency.pdf"),trnacodonabun_p,width=20,height=20)

for (phylum in unique(trna_codon_freq$PHYLUM)) {
  i_trnacodonabun_p <- ggplot(trna_codon_freq %>% filter(PHYLUM==phylum) %>% filter(! is.na(frequency)),
                          aes(color=SUBPHYLUM,fill=SUBPHYLUM)) + 
  geom_point(aes(x=relative_trna_abundance,y=frequency)) + theme_cowplot(12)  + facet_wrap(~trna_codon,nrow=8,ncol=8) +
  scale_colour_brewer(palette = "Set1") + ylab("codon abundance") + xlab("relative tRNA abundance") + 
    ggtitle(sprintf("%s tRNA abundance plot",phylum))
  ggsave(file.path(statsplotdir,sprintf("tRNA_abundance_codon_frequency_%s.pdf",phylum)),i_trnacodonabun_p,width=20,height=20)
}

# closeup shop
dbDisconnect(con, shutdown = TRUE)
dev.off()

