#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser(description='An executible script to split a seurat object by patient and save as individual seurat objects.')
parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object", help="Provide the full path to the seurat object for ehich you will run cellphoneDB on.")
# parser$add_argument("-o", "--output_directory", type="character", dest="output_directory", help="Provide the full path to directory in which to store the putput of this script.")
parser$add_argument("-id", "--sample_id_meta", type="character", dest="sample_id_meta", help="Provide the metadata column of the seurat object that contains the sampleID information.")

args <- parser$parse_args()
seurat_object <- args$seurat_object
# output_directory <- args$output_directory
sample_id_meta <- args$sample_id_meta

# #! del these
# seurat_object <- '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/adj_normal_subset/for_manu/draft3/highlevel/highevel_with_luad_matched_labels.Rda'
# sample_id_meta <- 'SampleID'
# output_directory <- '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/adj_normal_subset/for_manu/draft3/highlevel/'
# #!


library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)

srt <- get(load(file = seurat_object))

# split by sample ID
all.bypatient <- SplitObject(srt, split.by = sample_id_meta)

cpdb.output.path <- "cellphoneDB/patient_samples/"

if (!file.exists(file.path(cpdb.output.path))){
  dir.create(file.path(cpdb.output.path), showWarnings = FALSE, recursive = TRUE)
}

for (i in 1:length(all.bypatient)){
  saveRDS(all.bypatient[[i]], file = paste0(cpdb.output.path, names(all.bypatient[i]), ".rds"))
}


#* Maybe i could just output the folder where all of them are saved and then use somthing like this:
# my_pipeline( channel.from('/some/data') )