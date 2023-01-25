#######################################
###      Cell-type enrichment       ###
#######################################
### Code for cell-type enrichment was written by Suzanne Miedema
### Note that EWCE requires installation of additional R packages
# Create output directories
dir.create('figures/EWCE')
dir.create('output/EWCE')

library(readr)
library(readxl)
library(plyr)
library(data.table)
library(gplots)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)

##EWCE packages
library(EWCE)
library(cowplot)
library(limma)


####
#### input for EWCE
####

# # counts matrix = mat
# # cell type table = ct
# 
# # frontal cortical data
# mat.Lake <- read.table("lakeFrontal_norm_counts_all.txt", sep="\t", header=T, row.names = 1)
# ct.Lake <- read.table("lakeFrontal_celltypes.txt", sep="\t", header=T, row.names = 1)
# 
# 
# ####
# #### EWCE procedure - Specificity Matrix
# ####
# 
# ## Lake (Frontal)
# 
# # Make specificity matrix
# mat.Lake <- t(mat.Lake)
# annotLevel.Lake <- list(l1=ct.Lake$Celltype)
# ct.Lake_data <- generate.celltype.data(exp=mat.Lake, annotLevels = annotLevel.Lake,
#                                        groupName = "Lake", no_cores=1)
load("data/CellTypeData_Lake.rda")

# specificity.Lake <- rownames_to_column(data.frame(ctd[[1]][["specificity"]]), "Genes")
# write_csv(specificity.Lake, "EWCE_Lake_specificitymatrix.csv")


####
#### EWCE procedure - Statistics NDC vs TAU FRONTAL CORTEX
####

# Perform EWCE on input gene/prot list + bg list
bg.prodia.tau <- read_excel("data/ALLprot_preppedforEWCE.xlsx")
hits.prodia.tau <- read_excel("data/SIGNTAUprot_preppedforEWCE.xlsx")
hits.prodia.tau.UP <- filter(hits.prodia.tau, hits.prodia.tau$log2FC_TAU_control > 0)
hits.prodia.tau.DOWN <- filter(hits.prodia.tau, hits.prodia.tau$log2FC_TAU_control < 0)

# EWCE ON ALL SIGN GENES (UP&DOWN samen)
# (guidelines say >10.000 reps for publishing)
# LET OP: zorg dat de juiste ctd data in je environment staat (in dit geval Lake data)!
res.prodia.tau <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau$HGNC_symbol, annotLevel=1,
                                      bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                      sctSpecies="human", reps=20000) #, 
#res.prodia.tau <- bootstrap.enrichment.test(sct_data=ctd, human.hits=hits.prodia.tau$HGNC_symbol, #annotLevel=1,
#                                            human.bg=bg.prodia.tau$HGNC_symbol, geneSizeControl=TRUE, #genelistSpecies="human",
#                                            reps=20000) #, sctSpecies="human"

# EWCE ON ALL SIGN GENES (UP vs DOWN)
# (guidelines say >10.000 reps for publishing)
# Did this myself, in 2 steps, instead of using ewce_expression_data, as this gave me errors
res.prodia.tau.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.UP$HGNC_symbol, annotLevel=1,
                                       bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                       sctSpecies="human", reps=20000)

res.prodia.tau.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.DOWN$HGNC_symbol, annotLevel=1,
                                         bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                         sctSpecies="human", reps=20000)

# Save all results
write_csv(res.prodia.tau$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_SIGNPROTALL_BasedonLake.csv")
write_csv(res.prodia.tau.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_SIGNPROTUP_BasedonLake.csv")
write_csv(res.prodia.tau.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_SIGNPROTDOWN_BasedonLake.csv")

# Figure UP vs DOWN graph
merged_results.prodia.tau <- rbind(data.frame(res.prodia.tau.up$results,list="up"),
                                   data.frame(res.prodia.tau.down$results,list="down"))

pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau, mtc_method="BH")$plain
dev.off()



####
#### EWCE procedure - WGCNA MODULES - PER MODULE
####

# Perform EWCE on input gene/prot list + bg list
bg.prodia.tau <- read_excel("data/ALLprot_preppedforEWCE.xlsx")
hits.prodia.tau.ALLCLUST <- read_excel("data/211201_modules_prepped.xlsx")
hits.prodia.tau.blue <- filter(hits.prodia.tau.ALLCLUST, module=="blue")
hits.prodia.tau.blue.UP <- filter(hits.prodia.tau.blue, hits.prodia.tau.blue$log2FC_TAU_control > 0)
hits.prodia.tau.blue.DOWN <- filter(hits.prodia.tau.blue, hits.prodia.tau.blue$log2FC_TAU_control < 0)
hits.prodia.tau.brown <- filter(hits.prodia.tau.ALLCLUST, module=="brown")
hits.prodia.tau.brown.UP <- filter(hits.prodia.tau.brown, hits.prodia.tau.brown$log2FC_TAU_control > 0)
hits.prodia.tau.brown.DOWN <- filter(hits.prodia.tau.brown, hits.prodia.tau.brown$log2FC_TAU_control < 0)
hits.prodia.tau.yellow <- filter(hits.prodia.tau.ALLCLUST, module=="yellow")
hits.prodia.tau.yellow.UP <- filter(hits.prodia.tau.yellow, hits.prodia.tau.yellow$log2FC_TAU_control > 0)
hits.prodia.tau.yellow.DOWN <- filter(hits.prodia.tau.yellow, hits.prodia.tau.yellow$log2FC_TAU_control < 0)
hits.prodia.tau.pink <- filter(hits.prodia.tau.ALLCLUST, module=="pink")
hits.prodia.tau.pink.UP <- filter(hits.prodia.tau.pink, hits.prodia.tau.pink$log2FC_TAU_control > 0)
hits.prodia.tau.pink.DOWN <- filter(hits.prodia.tau.pink, hits.prodia.tau.pink$log2FC_TAU_control < 0)
hits.prodia.tau.black <- filter(hits.prodia.tau.ALLCLUST, module=="black")
hits.prodia.tau.black.UP <- filter(hits.prodia.tau.black, hits.prodia.tau.black$log2FC_TAU_control > 0)
hits.prodia.tau.black.DOWN <- filter(hits.prodia.tau.black, hits.prodia.tau.black$log2FC_TAU_control < 0)
hits.prodia.tau.turquoise <- filter(hits.prodia.tau.ALLCLUST, module=="turquoise")
hits.prodia.tau.turquoise.UP <- filter(hits.prodia.tau.turquoise, hits.prodia.tau.turquoise$log2FC_TAU_control > 0)
hits.prodia.tau.turquoise.DOWN <- filter(hits.prodia.tau.turquoise, hits.prodia.tau.turquoise$log2FC_TAU_control < 0)
hits.prodia.tau.lightcyan <- filter(hits.prodia.tau.ALLCLUST, module=="lightcyan")
hits.prodia.tau.lightcyan.UP <- filter(hits.prodia.tau.lightcyan, hits.prodia.tau.lightcyan$log2FC_TAU_control > 0)
hits.prodia.tau.lightcyan.DOWN <- filter(hits.prodia.tau.lightcyan, hits.prodia.tau.lightcyan$log2FC_TAU_control < 0)
hits.prodia.tau.green <- filter(hits.prodia.tau.ALLCLUST, module=="green")
hits.prodia.tau.green.UP <- filter(hits.prodia.tau.green, hits.prodia.tau.green$log2FC_TAU_control > 0)
hits.prodia.tau.green.DOWN <- filter(hits.prodia.tau.green, hits.prodia.tau.green$log2FC_TAU_control < 0)
hits.prodia.tau.red <- filter(hits.prodia.tau.ALLCLUST, module=="red")
hits.prodia.tau.red.UP <- filter(hits.prodia.tau.red, hits.prodia.tau.red$log2FC_TAU_control > 0)
hits.prodia.tau.red.DOWN <- filter(hits.prodia.tau.red, hits.prodia.tau.red$log2FC_TAU_control < 0)
hits.prodia.tau.greenyellow <- filter(hits.prodia.tau.ALLCLUST, module=="greenyellow")
hits.prodia.tau.greenyellow.UP <- filter(hits.prodia.tau.greenyellow, hits.prodia.tau.greenyellow$log2FC_TAU_control > 0)
hits.prodia.tau.greenyellow.DOWN <- filter(hits.prodia.tau.greenyellow, hits.prodia.tau.greenyellow$log2FC_TAU_control < 0)
hits.prodia.tau.tan <- filter(hits.prodia.tau.ALLCLUST, module=="tan")
hits.prodia.tau.tan.UP <- filter(hits.prodia.tau.tan, hits.prodia.tau.tan$log2FC_TAU_control > 0)
hits.prodia.tau.tan.DOWN <- filter(hits.prodia.tau.tan, hits.prodia.tau.tan$log2FC_TAU_control < 0)
hits.prodia.tau.purple <- filter(hits.prodia.tau.ALLCLUST, module=="purple")
hits.prodia.tau.purple.UP <- filter(hits.prodia.tau.purple, hits.prodia.tau.purple$log2FC_TAU_control > 0)
hits.prodia.tau.purple.DOWN <- filter(hits.prodia.tau.purple, hits.prodia.tau.purple$log2FC_TAU_control < 0)
hits.prodia.tau.magenta <- filter(hits.prodia.tau.ALLCLUST, module=="magenta")
hits.prodia.tau.magenta.UP <- filter(hits.prodia.tau.magenta, hits.prodia.tau.magenta$log2FC_TAU_control > 0)
hits.prodia.tau.magenta.DOWN <- filter(hits.prodia.tau.magenta, hits.prodia.tau.magenta$log2FC_TAU_control < 0)
hits.prodia.tau.grey60 <- filter(hits.prodia.tau.ALLCLUST, module=="grey60")
hits.prodia.tau.grey60.UP <- filter(hits.prodia.tau.grey60, hits.prodia.tau.grey60$log2FC_TAU_control > 0)
hits.prodia.tau.grey60.DOWN <- filter(hits.prodia.tau.grey60, hits.prodia.tau.grey60$log2FC_TAU_control < 0)
hits.prodia.tau.cyan <- filter(hits.prodia.tau.ALLCLUST, module=="cyan")
hits.prodia.tau.cyan.UP <- filter(hits.prodia.tau.cyan, hits.prodia.tau.cyan$log2FC_TAU_control > 0)
hits.prodia.tau.cyan.DOWN <- filter(hits.prodia.tau.cyan, hits.prodia.tau.cyan$log2FC_TAU_control < 0)
hits.prodia.tau.grey <- filter(hits.prodia.tau.ALLCLUST, module=="grey")
hits.prodia.tau.grey.UP <- filter(hits.prodia.tau.grey, hits.prodia.tau.grey$log2FC_TAU_control > 0)
hits.prodia.tau.grey.DOWN <- filter(hits.prodia.tau.grey, hits.prodia.tau.grey$log2FC_TAU_control < 0)
hits.prodia.tau.salmon <- filter(hits.prodia.tau.ALLCLUST, module=="salmon")
hits.prodia.tau.salmon.UP <- filter(hits.prodia.tau.salmon, hits.prodia.tau.salmon$log2FC_TAU_control > 0)
hits.prodia.tau.salmon.DOWN <- filter(hits.prodia.tau.salmon, hits.prodia.tau.salmon$log2FC_TAU_control < 0)
hits.prodia.tau.midnightblue <- filter(hits.prodia.tau.ALLCLUST, module=="midnightblue")
hits.prodia.tau.midnightblue.UP <- filter(hits.prodia.tau.midnightblue, hits.prodia.tau.midnightblue$log2FC_TAU_control > 0)
hits.prodia.tau.midnightblue.DOWN <- filter(hits.prodia.tau.midnightblue, hits.prodia.tau.midnightblue$log2FC_TAU_control < 0)


# EWCE ON GENES (UP vs DOWN) - PER MODULE
# (guidelines say >10.000 reps for publishing)
# Did this myself, in 2 steps, instead of using ewce_expression_data, as this gave me errors
res.prodia.blue.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.blue.UP$HGNC.symbol, annotLevel=1,
                                               bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                               sctSpecies="human", reps=20000)
res.prodia.blue.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.blue.DOWN$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.brown.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.brown.UP$HGNC.symbol, annotLevel=1,
                                                bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                sctSpecies="human", reps=20000)
res.prodia.brown.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.brown.DOWN$HGNC.symbol, annotLevel=1,
                                                  bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                  sctSpecies="human", reps=20000)
res.prodia.yellow.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.yellow.UP$HGNC.symbol, annotLevel=1,
                                                bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                sctSpecies="human", reps=20000)
res.prodia.yellow.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.yellow.DOWN$HGNC.symbol, annotLevel=1,
                                                  bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                  sctSpecies="human", reps=20000)
res.prodia.pink.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.pink.UP$HGNC.symbol, annotLevel=1,
                                                bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                sctSpecies="human", reps=20000)
res.prodia.pink.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.pink.DOWN$HGNC.symbol, annotLevel=1,
                                                  bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                  sctSpecies="human", reps=20000)
res.prodia.black.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.black.UP$HGNC.symbol, annotLevel=1,
                                                bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                sctSpecies="human", reps=20000)
res.prodia.black.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.black.DOWN$HGNC.symbol, annotLevel=1,
                                                  bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                  sctSpecies="human", reps=20000)
res.prodia.turquoise.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.turquoise.UP$HGNC.symbol, annotLevel=1,
                                                bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                sctSpecies="human", reps=20000)
res.prodia.turquoise.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.turquoise.DOWN$HGNC.symbol, annotLevel=1,
                                                  bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                  sctSpecies="human", reps=20000)
res.prodia.lightcyan.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.lightcyan.UP$HGNC.symbol, annotLevel=1,
                                                bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                sctSpecies="human", reps=20000)
res.prodia.lightcyan.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.lightcyan.DOWN$HGNC.symbol, annotLevel=1,
                                                  bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                  sctSpecies="human", reps=20000)
res.prodia.green.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.green.UP$HGNC.symbol, annotLevel=1,
                                                bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                sctSpecies="human", reps=20000)
res.prodia.green.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.green.DOWN$HGNC.symbol, annotLevel=1,
                                                  bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                  sctSpecies="human", reps=20000)
res.prodia.red.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.red.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.red.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.red.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
res.prodia.greenyellow.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.greenyellow.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.greenyellow.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.greenyellow.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
res.prodia.tan.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.tan.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.tan.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.tan.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
res.prodia.purple.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.purple.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.purple.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.purple.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
res.prodia.magenta.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.magenta.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.magenta.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.magenta.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
res.prodia.grey60.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.grey60.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.grey60.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.grey60.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
res.prodia.cyan.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.cyan.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.cyan.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.cyan.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
res.prodia.grey.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.grey.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.grey.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.grey.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
res.prodia.salmon.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.salmon.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.salmon.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.salmon.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
res.prodia.midnightblue.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.midnightblue.UP$HGNC.symbol, annotLevel=1,
                                                 bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                 sctSpecies="human", reps=20000)
res.prodia.midnightblue.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.midnightblue.DOWN$HGNC.symbol, annotLevel=1,
                                                   bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                   sctSpecies="human", reps=20000)
# Save all results
write_csv(res.prodia.blue.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-blue_UP_BasedonLake.csv")
write_csv(res.prodia.blue.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-blue_DOWN_BasedonLake.csv")
write_csv(res.prodia.brown.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-brown_UP_BasedonLake.csv")
write_csv(res.prodia.brown.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-brown_DOWN_BasedonLake.csv")
write_csv(res.prodia.yellow.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-yellow_UP_BasedonLake.csv")
write_csv(res.prodia.yellow.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-yellow_DOWN_BasedonLake.csv")
write_csv(res.prodia.pink.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-pink_UP_BasedonLake.csv")
write_csv(res.prodia.pink.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-pink_DOWN_BasedonLake.csv")
write_csv(res.prodia.black.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-black_UP_BasedonLake.csv")
write_csv(res.prodia.black.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-black_DOWN_BasedonLake.csv")
write_csv(res.prodia.turquoise.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-turquoise_UP_BasedonLake.csv")
write_csv(res.prodia.turquoise.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-turquoise_DOWN_BasedonLake.csv")
write_csv(res.prodia.lightcyan.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-lightcyan_UP_BasedonLake.csv")
write_csv(res.prodia.lightcyan.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-lightcyan_DOWN_BasedonLake.csv")
write_csv(res.prodia.green.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-green_UP_BasedonLake.csv")
write_csv(res.prodia.green.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-green_DOWN_BasedonLake.csv")
write_csv(res.prodia.red.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-red_UP_BasedonLake.csv")
write_csv(res.prodia.red.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-red_DOWN_BasedonLake.csv")
write_csv(res.prodia.greenyellow.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-greenyellow_UP_BasedonLake.csv")
write_csv(res.prodia.greenyellow.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-greenyellow_DOWN_BasedonLake.csv")
write_csv(res.prodia.tan.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-tan_UP_BasedonLake.csv")
write_csv(res.prodia.tan.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-tan_DOWN_BasedonLake.csv")
write_csv(res.prodia.purple.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-purple_UP_BasedonLake.csv")
write_csv(res.prodia.purple.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-purple_DOWN_BasedonLake.csv")
write_csv(res.prodia.magenta.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-magenta_UP_BasedonLake.csv")
write_csv(res.prodia.magenta.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-magenta_DOWN_BasedonLake.csv")
write_csv(res.prodia.grey60.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-grey60_UP_BasedonLake.csv")
write_csv(res.prodia.grey60.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-grey60_DOWN_BasedonLake.csv")
write_csv(res.prodia.cyan.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-cyan_UP_BasedonLake.csv")
write_csv(res.prodia.cyan.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-cyan_DOWN_BasedonLake.csv")
write_csv(res.prodia.grey.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-grey_UP_BasedonLake.csv")
write_csv(res.prodia.grey.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-grey_DOWN_BasedonLake.csv")
write_csv(res.prodia.salmon.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-salmon_UP_BasedonLake.csv")
write_csv(res.prodia.salmon.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-salmon_DOWN_BasedonLake.csv")
write_csv(res.prodia.midnightblue.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-midnightblue_UP_BasedonLake.csv")
write_csv(res.prodia.midnightblue.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-midnightblue_DOWN_BasedonLake.csv")

# Figures UP vs DOWN graph
merged_results.prodia.tau.blue <- rbind(data.frame(res.prodia.blue.up$results,list="up"),
                                   data.frame(res.prodia.blue.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-blue_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.blue, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.brown <- rbind(data.frame(res.prodia.brown.up$results,list="up"),
                                        data.frame(res.prodia.brown.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-brown_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.brown, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.yellow <- rbind(data.frame(res.prodia.yellow.up$results,list="up"),
                                        data.frame(res.prodia.yellow.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-yellow_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.yellow, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.pink <- rbind(data.frame(res.prodia.pink.up$results,list="up"),
                                        data.frame(res.prodia.pink.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-pink_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.pink, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.black <- rbind(data.frame(res.prodia.black.up$results,list="up"),
                                        data.frame(res.prodia.black.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-black_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.black, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.turquoise <- rbind(data.frame(res.prodia.turquoise.up$results,list="up"),
                                        data.frame(res.prodia.turquoise.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-turquoise_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.turquoise, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.lightcyan <- rbind(data.frame(res.prodia.lightcyan.up$results,list="up"),
                                        data.frame(res.prodia.lightcyan.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-lightcyan_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.lightcyan, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.green <- rbind(data.frame(res.prodia.green.up$results,list="up"),
                                        data.frame(res.prodia.green.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-green_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.green, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.red <- rbind(data.frame(res.prodia.red.up$results,list="up"),
                                        data.frame(res.prodia.red.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-red_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.red, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.greenyellow <- rbind(data.frame(res.prodia.greenyellow.up$results,list="up"),
                                        data.frame(res.prodia.greenyellow.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-greenyellow_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.greenyellow, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.tan <- rbind(data.frame(res.prodia.tan.up$results,list="up"),
                                        data.frame(res.prodia.tan.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-tan_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.tan, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.purple <- rbind(data.frame(res.prodia.purple.up$results,list="up"),
                                        data.frame(res.prodia.purple.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-purple_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.purple, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.magenta <- rbind(data.frame(res.prodia.magenta.up$results,list="up"),
                                        data.frame(res.prodia.magenta.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-magenta_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.magenta, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.grey60 <- rbind(data.frame(res.prodia.grey60.up$results,list="up"),
                                        data.frame(res.prodia.grey60.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-grey60_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.grey60, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.cyan <- rbind(data.frame(res.prodia.cyan.up$results,list="up"),
                                        data.frame(res.prodia.cyan.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-cyan_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.cyan, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.grey <- rbind(data.frame(res.prodia.grey.up$results,list="up"),
                                        data.frame(res.prodia.grey.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-grey_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.grey, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.salmon <- rbind(data.frame(res.prodia.salmon.up$results,list="up"),
                                                 data.frame(res.prodia.salmon.down$results,list="down"))

pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-salmon_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.salmon, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.midnightblue <- rbind(data.frame(res.prodia.midnightblue.up$results,list="up"),
                                          data.frame(res.prodia.midnightblue.down$results,list="down"))

pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-midnightblue_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.midnightblue, mtc_method="BH")$plain
dev.off()


####
#### EWCE procedure - WGCNA MODULES - FRONTAL CORTEX VALIDATED MODULES VS NON-VALIDATED MODULES
####

# The modules that were validated in the RiMOD frontal cortex FTD-MAPT data set were: 
# Blue, Brown, Cyan, Greenyellow, Red, Purple, Turquoise, Yellow
# The modules that were NOT validated in the RiMOD frontal cortex FTD-MAPT data set were: 
# Black, Green, Grey, Grey60, Magenta, Midnightblue, Lightcyan, Pink, Salmon, Tan

# Perform EWCE on input gene/prot list + bg list
bg.prodia.tau <- read_excel("data/ALLprot_preppedforEWCE.xlsx")
hits.prodia.tau.ALLCLUST <- read_excel("data/211201_modules_prepped.xlsx")
hits.prodia.tau.FRONT.VALID <- filter(hits.prodia.tau.ALLCLUST, module %in%
                                        c("turquoise","blue","brown","yellow","red","purple","greenyellow","cyan"))
hits.prodia.tau.FRONT.VALID.UP <- filter(hits.prodia.tau.FRONT.VALID, hits.prodia.tau.FRONT.VALID$log2FC_TAU_control > 0)
hits.prodia.tau.FRONT.VALID.DOWN <- filter(hits.prodia.tau.FRONT.VALID, hits.prodia.tau.FRONT.VALID$log2FC_TAU_control < 0)
hits.prodia.tau.FRONT.NOTVALID <- filter(hits.prodia.tau.ALLCLUST, module %in%
                                           c("grey","black","pink","magenta","salmon","midnightblue","lightcyan","grey60","green","tan"))
hits.prodia.tau.FRONT.NOTVALID.UP <- filter(hits.prodia.tau.FRONT.NOTVALID, hits.prodia.tau.FRONT.NOTVALID$log2FC_TAU_control > 0)
hits.prodia.tau.FRONT.NOTVALID.DOWN <- filter(hits.prodia.tau.FRONT.NOTVALID, hits.prodia.tau.FRONT.NOTVALID$log2FC_TAU_control < 0)


# EWCE ON GENES (UP vs DOWN) - PER VALIDATED and NON-VALIDATED modules FRONTAL CORTEX
# (guidelines say >10.000 reps for publishing)
# Did this myself, in 2 steps, instead of using ewce_expression_data, as this gave me errors
res.prodia.tau.FRONT.VALID.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.FRONT.VALID.UP$HGNC.symbol, annotLevel=1,
                                                bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                sctSpecies="human", reps=20000)
res.prodia.tau.FRONT.VALID.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.FRONT.VALID.DOWN$HGNC.symbol, annotLevel=1,
                                                  bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                  sctSpecies="human", reps=20000)
res.prodia.tau.FRONT.NOTVALID.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.FRONT.NOTVALID.UP$HGNC.symbol, annotLevel=1,
                                                           bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                           sctSpecies="human", reps=20000)
res.prodia.tau.FRONT.NOTVALID.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.FRONT.NOTVALID.DOWN$HGNC.symbol, annotLevel=1,
                                                             bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                             sctSpecies="human", reps=20000)

# Save all results
write_csv(res.prodia.tau.FRONT.VALID.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-ValidatedinFrontCortRiMODFTD-MAPT_UP_BasedonLake.csv")
write_csv(res.prodia.tau.FRONT.VALID.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-ValidatedinFrontCortRiMODFTD-MAPT_DOWN_BasedonLake.csv")
write_csv(res.prodia.tau.FRONT.NOTVALID.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-NOTValidatedinFrontCortRiMODFTD-MAPT_UP_BasedonLake.csv")
write_csv(res.prodia.tau.FRONT.NOTVALID.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-NOTValidatedinFrontCortRiMODFTD-MAPT_DOWN_BasedonLake.csv")

# Figures UP vs DOWN graph
merged_results.prodia.tau.FRONT.VALID <- rbind(data.frame(res.prodia.tau.FRONT.VALID.up$results,list="up"),
                                        data.frame(res.prodia.tau.FRONT.VALID.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-ValidatedinFrontCortRiMODFTD-MAPT_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.FRONT.VALID, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.FRONT.NOTVALID <- rbind(data.frame(res.prodia.tau.FRONT.NOTVALID.up$results,list="up"),
                                         data.frame(res.prodia.tau.FRONT.NOTVALID.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-NOTValidatedinFrontCortRiMODFTD-MAPT_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.FRONT.NOTVALID, mtc_method="BH")$plain
dev.off()


####
#### EWCE procedure - WGCNA MODULES - TEMPORAL CORTEX VALIDATED MODULES VS NON-VALIDATED MODULES
####

# The modules that were validated in the RiMOD temporal cortex FTD-MAPT data set were: 
# Black, Blue, Brown, Cyan, Green, Magenta, Purple, Red, Turquoise, Yellow 
# The modules that were NOT validated in the RiMOD temporal cortex FTD-MAPT data set were: 
# Greenyellow, Grey, Grey60, Lightcyan, Midnightblue, Pink, Salmon, Tan

# Perform EWCE on input gene/prot list + bg list
bg.prodia.tau <- read_excel("data/ALLprot_preppedforEWCE.xlsx")
hits.prodia.tau.ALLCLUST <- read_excel("data/211201_modules_prepped.xlsx")
hits.prodia.tau.TEMP.VALID <- filter(hits.prodia.tau.ALLCLUST, module %in%
                                        c("turquoise","blue","brown","yellow","green","red","black","purple","cyan","magenta"))
hits.prodia.tau.TEMP.VALID.UP <- filter(hits.prodia.tau.TEMP.VALID, hits.prodia.tau.TEMP.VALID$log2FC_TAU_control > 0)
hits.prodia.tau.TEMP.VALID.DOWN <- filter(hits.prodia.tau.TEMP.VALID, hits.prodia.tau.TEMP.VALID$log2FC_TAU_control < 0)
hits.prodia.tau.TEMP.NOTVALID <- filter(hits.prodia.tau.ALLCLUST, module %in%
                                           c("grey","greenyellow","tan","salmon","midnightblue","lightcyan","grey60","pink"))
hits.prodia.tau.TEMP.NOTVALID.UP <- filter(hits.prodia.tau.TEMP.NOTVALID, hits.prodia.tau.TEMP.NOTVALID$log2FC_TAU_control > 0)
hits.prodia.tau.TEMP.NOTVALID.DOWN <- filter(hits.prodia.tau.TEMP.NOTVALID, hits.prodia.tau.TEMP.NOTVALID$log2FC_TAU_control < 0)


# EWCE ON GENES (UP vs DOWN) - PER VALIDATED and NON-VALIDATED modules TEMPORAL CORTEX
# (guidelines say >10.000 reps for publishing)
# Did this myself, in 2 steps, instead of using ewce_expression_data, as this gave me errors
res.prodia.tau.TEMP.VALID.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.TEMP.VALID.UP$HGNC.symbol, annotLevel=1,
                                                           bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                           sctSpecies="human", reps=20000)
res.prodia.tau.TEMP.VALID.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.TEMP.VALID.DOWN$HGNC.symbol, annotLevel=1,
                                                             bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                             sctSpecies="human", reps=20000)
res.prodia.tau.TEMP.NOTVALID.up <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.TEMP.NOTVALID.UP$HGNC.symbol, annotLevel=1,
                                                              bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                              sctSpecies="human", reps=20000)
res.prodia.tau.TEMP.NOTVALID.down <- bootstrap.enrichment.test(sct_data=ctd, hits=hits.prodia.tau.TEMP.NOTVALID.DOWN$HGNC.symbol, annotLevel=1,
                                                                bg=bg.prodia.tau$HGNC_symbol, genelistSpecies="human",
                                                                sctSpecies="human", reps=20000)

# Save all results
write_csv(res.prodia.tau.TEMP.VALID.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-ValidatedinTempCortRiMODFTD-MAPT_UP_BasedonLake.csv")
write_csv(res.prodia.tau.TEMP.VALID.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-ValidatedinTempCortRiMODFTD-MAPT_DOWN_BasedonLake.csv")
write_csv(res.prodia.tau.TEMP.NOTVALID.up$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-NOTValidatedinTempCortRiMODFTD-MAPT_UP_BasedonLake.csv")
write_csv(res.prodia.tau.TEMP.NOTVALID.down$results, "output/EWCE/EWCE_frontal_PRODIA_TAU_M-NOTValidatedinTempCortRiMODFTD-MAPT_DOWN_BasedonLake.csv")

# Figures UP vs DOWN graph
merged_results.prodia.tau.TEMP.VALID <- rbind(data.frame(res.prodia.tau.TEMP.VALID.up$results,list="up"),
                                         data.frame(res.prodia.tau.TEMP.VALID.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-ValidatedinTempCortRiMODFTD-MAPT_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.TEMP.VALID, mtc_method="BH")$plain
dev.off()

merged_results.prodia.tau.TEMP.NOTVALID <- rbind(data.frame(res.prodia.tau.TEMP.NOTVALID.up$results,list="up"),
                                            data.frame(res.prodia.tau.TEMP.NOTVALID.down$results,list="down"))
pdf(file = "figures/EWCE/EWCE_frontal_PRODIA_TAU_M-NOTValidatedinTempCortRiMODFTD-MAPT_UP&DOWN_BasedonLake.pdf", height = 7.5, width = 9.5)
ewce.plot(total_res=merged_results.prodia.tau.TEMP.NOTVALID, mtc_method="BH")$plain
dev.off()
