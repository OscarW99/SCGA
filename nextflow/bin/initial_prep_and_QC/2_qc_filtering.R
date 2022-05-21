#!/usr/bin/env Rscript#!/usr/bin/env Rscript


library(argparse)
library(Seurat)
library(rjson)
library(ggplot2)
library(sctransform)

parser <- ArgumentParser(description='An executible R script to filter cells in a seurta object based on provided thresholds for nCounts, nFeatures and percent.mt genes.')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object_path", help="Provide the full path to the seuart object .Rda file.")
parser$add_argument("-pf", "--parameter_file", type="character", dest="parameter_file", help="Provide the path to a json file that contains threshold values for nCounts, nFeatures and percent.mt genes.")

args <- parser$parse_args()
seurat_object_path <- args$seurat_object_path
parameter_file <- args$parameter_file

srt <- get(load(file=seurat_object_path))
parameters <- fromJSON(file=parameter_file)

# # # Example JSON file format
# {   
#     "QC_thresholds":{
#         "nFeatures_threshold":400,
#         "nCounts_threshold":1000,
#         "percent_mt_threshold":25,
#         "cluster_resolution":0.5
#         },
#     "Gene_sets":{
#         "Myeloid":['LYZ', 'ABC', 'WTF']
#         'Endothelial':['ect', 'lol', 'omg']
#     }
    

# }

cells_before_filter <- length(colnames(srt))

# QC filtering
srt <- subset(srt, subset = nFeature_RNA > as.integer(parameters$nFeatures_threshold) & nCount_RNA > as.integer(parameters$nCounts_threshold) & percent.mt < as.integer(parameters$percent_mt_threshold))

cells_after_filter <- length(colnames(srt))

# Regressing out the difference between the G2M and S phase scores. This means that signals separating non-cycling cells and cycling cells will be maintained, but differences in cell cycle phase amongst proliferating cells (which are often uninteresting), will be regressed out of the data
# marrow$CC.Difference <- marrow$S.Score - marrow$G2M.Score
# marrow <- ScaleData(marrow, vars.to.regress = "CC.Difference", features = rownames(marrow))
# marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)


options(future.globals.maxSize = 4000 * 1024^2)
srt <- SCTransform(srt, vars.to.regress = "percent.mt", verbose = FALSE)
srt <- RunPCA(srt, verbose = FALSE)
srt <- RunUMAP(srt, dims = 1:30, verbose = FALSE)
srt <- FindNeighbors(srt, dims = 1:30, verbose = FALSE)
srt <- FindClusters(srt, verbose = FALSE, resolution = as.integer(parameters$cluster_resolution)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
srt <- CellCycleScoring(srt, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)


# Produce pots at this stage?
# 1) UMAPs for clusters, cell cycle score, %mt, nCount, nFeatures
# 2) Top markers for each cluster in the form of a heatmap

# 3) I will want to be able to play around with different plots etc. I think saving the srt object here and then 
#  The next script will be where I assign celltpye names to the cells


#*----------------------------------------------------------------------------------------------------------
# I should do somthing like those people at mt-sinai whos nextflow pipeline produces pdf slides with each figure explaining what they are. I could also then include what the statistics show ect.

#! How will I Output the number of cells filtered? As a graph, A table? both. Also is there an R package to create a publication quality figure. It would be good to auto-generate figures and figure legnends.

#? Based on the tissue type I could have known markers plotted with markers in this tsv file... https://panglaodb.se/markers.html?cell_type=%27choose%27
#? I may have to re-download this file periodically for updated versions? Although it hasn't been updated for 2 years so perhaps I wont need to do this.
#? I can also add my own markers to this file. So at least for the lung I can add them and then perhaps manually adding ones per tissue based on the literatire for specific clients. Eventually building up a marker database perhaps for mouse and human.