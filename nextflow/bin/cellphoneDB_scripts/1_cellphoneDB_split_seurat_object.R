#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser(description='An executible script to split a seurat object by patient and save as individual seurat objects.')
parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object", help="Provide the full path to the seurat object for ehich you will run cellphoneDB on.")
parser$add_argument("-id", "--sample_id_meta", type="character", dest="sample_id_meta", help="Provide the metadata column of the seurat object that contains the sampleID information.")
parser$add_argument("-l", "--celltype_label_meta", type="character", dest="celltype_label_meta", help="Provide the metadata column of the seurat object that contains the celltype/compartment information. This is the levels at which we will compare cell interactions.")

args <- parser$parse_args()
seurat_object <- args$seurat_object
sample_id_meta <- args$sample_id_meta
celltype_label_meta <- args$celltype_label_meta


# #! del these
# seurat_object <- '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/adj_normal_subset/for_manu/draft3/highlevel/highevel_with_luad_matched_labels.Rda'
# sample_id_meta <- 'sampleID'
# celltype_label_meta <- 'celltype_highlevel'
# output_directory <- '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/adj_normal_subset/for_manu/draft3/highlevel/'
# #!

library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)
library(future.apply)

srt <- get(load(file = seurat_object))

# split by sampleID_meta
all.bypatient <- SplitObject(srt, split.by = sample_id_meta)

options(future.globals.maxSize = 48000 * 1024^2)

# todo - could I use parallel computing practices within R to speed this bit up?
# todo - rename cpdb_prep (1st) process
for (i in 1:length(all.bypatient)){
  srt <- all.bypatient[[i]]
  # Create input files for cpdb
  sampleid <- unique(srt@meta.data[,sample_id_meta])
  # dir.create(file.path(sampleid), showWarnings = FALSE)
  counts <- srt[["RNA"]]@counts
  counts.norm <- future_apply(counts, 2, function(x) (x/sum(x))*10000) # this is recommended by the cellphonedb paper
  write.table(counts.norm, paste0(sampleid, "_counts.txt"), sep="\t", quote=F)
  meta <- srt@meta.data
  meta <- cbind(rownames(meta), meta[,celltype_label_meta, drop=F])
  meta[is.na(meta[,2]),2]<-'Unassigned'
  colnames(meta)<-c('cell','annotation')
  write.table(meta, paste0(sampleid, "_meta.txt"), sep="\t", quote=F, row.names=F)
}

# #* Testing ##
# counts.norm <- "hello world"
# meta <- "hello metaverse"
# sampleid <- "patient4236"
# write.table(counts.norm, paste0(sampleid, "_counts.txt"), sep="\t", quote=F)
# write.table(meta, paste0(sampleid, "_meta.txt"), sep="\t", quote=F, row.names=F)
# #*###########


# # The output of the process that hold this script will be the counts and txt files. 
# path "${x}_meta_highlvl.txt" into cpdb_file_prep
# path "${x}_counts_highlvl.txt" into cpdb_file_prep

# # The input to the next process (where cpdb is actually run) could look like this...
# cpdb_file_prep
#     .fromFilePairs('${workDir}/*_{meta,counts}_highlvl.txt')


# # https://www.nextflow.io/docs/latest/process.html#dynamic-output-file-names