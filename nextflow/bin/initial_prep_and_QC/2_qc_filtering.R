#!/usr/bin/env Rscript#!/usr/bin/env Rscript


library(argparse)
library(Seurat)
library(rjson)

parser <- ArgumentParser(description='An executible R script to filter cells in a seurta object based on provided thresholds for nCounts, nFeatures and percent.mt genes.')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object_path", help="Provide the full path to the seuart object .Rda file.")
parser$add_argument("-pf", "--parameter_file", type="character", dest="parameter_file", help="Provide the path to a json file that contains threshold values for nCounts, nFeatures and percent.mt genes.")

args <- parser$parse_args()
seurat_object_path <- args$seurat_object_path
parameter_file <- args$parameter_file

srt <- get(load(file=seurat_object_path))
parameters <- fromJSON(file=parameter_file)

# # JSON file format
# { 
#     "nFeatures_threshold":400,
#     "nCounts_threshold":1000,
#     "percent_mt_threshold":25
# }

# QC filtering
srt <- subset(srt, subset = nFeature_RNA > as.integer(parameters$nFeatures_threshold) & nCount_RNA > as.integer(parameters$nCounts_threshold) & percent.mt < as.integer(parameters$percent_mt_threshold))



s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
srt <- CellCycleScoring(srt, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)