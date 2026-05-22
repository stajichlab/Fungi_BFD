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
outdir="reports"
dir.create(file.path(outdir),showWarnings = FALSE)
pfamoutdir = file.path(outdir,"pfam")
dir.create(pfamoutdir, showWarnings = FALSE)

pdf(file.path(pfamoutdir,"pfam_sum.pdf"))
# to use a database file already created by 
con <- dbConnect(duckdb(), dbdir=file.path(DBDIR,DBNAME), read_only = TRUE)
pfamstat_sql ="
SELECT sp.*, ct.pfam_id, ct.pfam_count as pfam_count
FROM 
species AS sp
JOIN ( SELECT species_prefix, pfam_id, COUNT(*) AS pfam_count
FROM pfam GROUP BY pfam_id, species_prefix) AS ct
ON ct.species_prefix = sp.LOCUSTAG
"

pfamstat_res <- dbGetQuery(con, pfamstat_sql)
head(pfamstat_res)

# sanity check
# pfamstat_res %>% filter(pfam_id == "AAA_30")
# sanity check
#pfamstat_res %>%
#  count(pfam_id, SPECIESIN) %>%
#  filter(n > 1)

pivot <- pfamstat_res %>% 
  select(c(pfam_id, PHYLUM,SPECIESIN, pfam_count)) %>% 
  pivot_wider(names_from = c(SPECIESIN,PHYLUM), values_from = pfam_count,values_fill = 0)

write_csv(pivot,file.path(outdir,"pfam","pfam_counts.csv"))

crn <- pfamstat_res %>% 
  filter(pfam_id == "Crinkler" & pfam_count > 0) %>% 
  select(c(PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIESIN,pfam_count)) %>%
  arrange(PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIESIN)
write_csv(crn,file.path(outdir,"pfam","Crinkler_counts.csv"))

protease <- pfamstat_res %>% 
    filter(pfam_id == "Peptidase_M36" & pfam_count > 0) %>% 
  select(c(PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIESIN,pfam_count)) %>%
  arrange(PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIESIN)
write_csv(protease,file.path(outdir,"pfam","Peptidase_M36_counts.csv"))

