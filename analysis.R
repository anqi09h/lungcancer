getwd()
list.files()

library(BiocManager)
library(remotes)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyverse)
library(maftools)
library(DESeq2)
library(EnhancedVolcano)
library(dplyr)
library(ggplot2)
library(ggrepel)

gdc_project = getGDCprojects()
View(gdc_project)
head(gdc_project)
head(gdc_project)[c("project_id", "name")]

getGDCprojects()$project_id
getProjectSummary("TCGA-LUAD")

snv_query = GDCquery(project = 'TCGA-LUAD', 
                     data.category = 'Simple Nucleotide Variation')
View(snv_query)
output_snv_query = getResults(snv_query)
View(output_snv_query)
snv_query = GDCquery(project = 'TCGA-LUAD', 
                     data.category = 'Simple Nucleotide Variation',
                      access = 'open',
                     data.type = 'Masked Somatic Mutation',
                     data.format = 'MAF')
View(snv_query)

output_snv_query = getResults(snv_query)
View(output_snv_query)

GDCdownload(snv_query)
maf <- GDCprepare(snv_query)

maf = read.maf(maf = maf)

maf <- readRDS("exports/TCGA_LUAD_maf.rds")

View(maf)

getSampleSummary(maf)

plotmafSummary(maf = maf,
               rmOutlier = TRUE,
               dashboard = TRUE)

oncoplot(maf = maf,
         top = 20)

titv(maf = maf,
     plot = TRUE,
     useSyn = TRUE)
