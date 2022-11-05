#!/usr/bin/env Rscript

# For this one we will start by assuming the data will be in the form of count data for each patient and the folders will be named after the sampleID.

# So using Bischoff as an example for now. We wiil need to see what the actual output from 10X ect is. We could design the pipeline to integrate seamlessly from these sequencing technologies.
# This script currently requires that the folder names are exactly as below.

## EXAMPLE DIRECOTRY STRUCTURE ##
# |-- patient1
# |   `-- filtered_feature_bc_matrix
# |       |-- barcodes.tsv.gz
# |       |-- features.tsv.gz
# |       `-- matrix.mtx.gz
# |-- patient2
# |   `-- filtered_feature_bc_matrix
# |       |-- barcodes.tsv.gz
# |       |-- features.tsv.gz
# |       `-- matrix.mtx.gz
# `-- patient3
#     `-- filtered_feature_bc_matrix
#         |-- barcodes.tsv.gz
#         |-- features.tsv.gz
#         `-- matrix.mtx.gz


## Arguments?
# The root folder for these files
library(argparse)

parser <- ArgumentParser(description='An executible R script to create a merged seurat object from the single cell sequencing data from many patients.')

parser$add_argument("-d", "--data_directory", type="character", dest="data_directory", help="Provide the full path to directory in which all patient count data is stored.")

args <- parser$parse_args()
data_directory <- args$data_directory



# Import libraries
library(Seurat)
library(ggplot2)

# Loop through each patient folder to create sparce matrix for each patient
# 1)barcodes.tsv.gz    2)features.tsv.gz    3) matrix.mtx.gz 

patient_folders <- list.dirs(data_directory, full.names = TRUE, recursive = FALSE)

seurat_object_holder <- list()

for (i in 1:length(patient_folders)) {
    patient <- patient_folders[i]
	sparce_matrix <- Read10X(data.dir=paste0(patient_folders, "/filtered_feature_bc_matrix"))
	seurat_obj <- CreateSeuratObject(sparce_matrix)
	sample_id <- basename(patient)
    
	seurat_obj$sampleID <- sample_id
	seurat_obj$barcode <- colnames(seurat_obj)
	seurat_obj$percent.mt <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

	seurat_object_holder[[i]] <- seurat_obj
	print(paste0('creating seurat object for patient: ', sample_id))
}

soh <- unlist(seurat_object_holder, use.names=FALSE)

# loop o create merge command
#* add an 'if length soh is less then 3 option'
string <- "merge(soh[[1]], y = c(soh[[2]]"
for (i in 3:length(soh)){
	addition <- paste0(', soh[[', i, ']]')
	string <- paste0(string, addition)
}
string <- paste0(string, '))')

srt <- eval(parse(text=string))



# https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html
# The UMI counts per cell should generally be above 500, that is the low end of what we expect. If UMI counts are between 500-1000 counts, it is usable but the cells probably should have been sequenced more deeply.
# srt <- subset(srt, subset = nFeature_RNA > 400 & nCount_RNA > 1000 & percent.mt < 25) # percent.mt = v3:25, v2:10, v1:10
# Doublets cannot be accurately removed using feature counts and UMI counts (/ther are better ways). We will do this later with another tool.

#* GOOD FIGURE, USE FOR WORK
# library(ggplot2)
# library(dplyr)

# srt@meta.data %>% 
# 	ggplot(aes(color=dataset, x=nCount_RNA, fill= dataset)) + 
# 	geom_density(alpha = 0.2) + 
# 	scale_x_log10() + 
# 	theme_classic() +
# 	ylab("Cell density") +
# 	geom_vline(xintercept = 500)
# # VlnPlot(srt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3 , pt.size = 0)
# ggsave('/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/test.png')

save(srt, file="srt.Rda")

# Plot counts
ggplot(srt@meta.data, aes(x=nCount_RNA)) + 
geom_density(alpha = 0.2) + 
scale_x_log10() + 
theme_classic() +
ylab("Cell density") +
geom_vline(xintercept = 1000, color='red')
ggsave("nCount.png")

# Plot features
ggplot(srt@meta.data, aes(x=nFeature_RNA)) + 
geom_density(alpha = 0.2) + 
scale_x_log10() + 
theme_classic() +
ylab("Cell density") +
geom_vline(xintercept = 400, color='red')
ggsave("nFeature.png")

ggplot(srt@meta.data, aes(x=percent.mt)) + 
geom_density(alpha = 0.2) + 
scale_x_log10() + 
theme_classic() +
ylab("Cell density") +
geom_vline(xintercept = 25, color='red')
ggsave("percent.mt.png")

# To load the data again
#srt <- get(load(file="/ahg/regevdata/projects/lungCancerBueno/Results/10x_bischoff_102621/srt.Rda"))


# Output: (20680 cells filtered)
# ... (below is before filtering)
# p018n p018t p019n p019t p023t p024t p027n p027t p028n p029n p030n p030t p031n
# 12183 14770  1547  1557  7811  6753  9045  9371  6615  5173  2986  4861  5674
# p031t p032n p032t p033n p033t p034n p034t
#  6134  5900 11942  5739  5202  5308  5165


#@ RUNNING THE COMMAND LINE SCRIPT ##
# cd '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/PRIV_GITHUB/SCGA/nextflow/bin'
# Rscript 1_create_seurat_object.R -d '/ahg/regevdata/projects/lungCancerBueno/Results/10x_bischoff_102621/data/'
