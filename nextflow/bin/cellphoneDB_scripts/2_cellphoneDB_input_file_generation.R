#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description='An executible script to prepare files to run cellphoneDB in a patient-wise manner.')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object", help="Provide the full path to the seurat object for ehich you will run cellphoneDB on.")
parser$add_argument("-o", "--output_directory", type="character", dest="output_directory", help="Provide the full path to directory in which to store the putput of this script.")
parser$add_argument("-id", "--sample_id_meta", type="character", dest="sample_id_meta", help="Provide the metadata column of the seurat object that contains the sampleID information.")
parser$add_argument("-l", "--celltype_label_meta", type="character", dest="celltype_label_meta", help="Provide the metadata column of the seurat object that contains the celltype/compartment information. This is the levels at which we will compare cell interactions.")

args <- parser$parse_args()
seurat_object <- args$seurat_object
output_directory <- args$output_directory
sample_id_meta <- args$sample_id_meta
celltype_label_meta <- args$celltype_label_meta

#! del these
# seurat_object <- '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/adj_normal_subset/for_manu/draft3/highlevel/highevel_with_luad_matched_labels.Rda'
# sample_id_meta <- 'SampleID'
#celltype_label_meta <- 'luad_label_match'
# output_directory <- '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/adj_normal_subset/for_manu/draft3/highlevel/'
#!


library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(future.apply)



# LoadRDS
srt <- get(load(file = seurat_object))

# Load var from last script
cpdb.output.path <- paste0(output_directory, "/cellphoneDB/cpdb_input_files/")

if (!file.exists(file.path(cpdb.output.path))){
    dir.create(file.path(cpdb.output.path), showWarnings = FALSE, recursive = TRUE)
}

options(future.globals.maxSize = 48000 * 1024^2)

sampleid <- unique(srt@meta.data[,sample_id_meta])
dir.create(file.path(cpdb.output.path, sampleid), showWarnings = FALSE)
counts <- srt[["RNA"]]@counts
counts.norm <- future_apply(counts, 2, function(x) (x/sum(x))*10000) # this is recommended by the cellphonedb paper
write.table(counts.norm, paste0(cpdb.output.path, sampleid, "/", sampleid, "_counts.txt"), sep="\t", quote=F)
meta <- srt@meta.data
meta <- cbind(rownames(meta), meta[,celltype_label_meta, drop=F])
meta[is.na(meta[,2]),2]<-'Unassigned'
colnames(meta)<-c('cell','annotation')
write.table(meta, paste0(cpdb.output.path, sampleid, "/", sampleid, "_meta.txt"), sep="\t", quote=F, row.names=F)
