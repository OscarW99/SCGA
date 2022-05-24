#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description='An executible script to re-cluster a seurat object based on a given resolution.')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object_path", help="Provide the full path to the seuart object .Rda file.")
parser$add_argument("-pf", "--parameter_file", type="character", dest="parameter_file", help="Provide the path to a json file that contains clustering resolution")

#! The batch correction argument should be a string of all metadata columns seperated by a space. Also surround in quotes
args <- parser$parse_args()
seurat_object_path <- args$seurat_object_path
parameter_file <- args$parameter_file

message('Loading Seurat Object')
srt <- get(load(file=seurat_object_path))
parameters <- fromJSON(file=parameter_file)

# What json import will look like
# {   
#     "QC_thresholds":{
#         "nFeatures_threshold":400,
#         "nCounts_threshold":1000,
#         "percent_mt_threshold":25,
#         "cluster_resolution":0.5
#         },
#     "Gene_sets":{
#         "Myeloid":["S100A8", "RETN", "FCN1"],
#         "Bplasmast":["IGHGP", "IGHG1", "IGHG3"]
#     },
#     "Clusters":{
#         "number_of_clusters":5
#         "new_cluster_names":c("Tcell", "Bcel", "Myeloid", "Epithelial", "Endothelial")
#     }
# }

new_labels <- parameters$Clusters$new_cluster_names
srt$compartment <- lapply(srt$seurat_clusters, function(x) new_labels[x+1])

DimPlot(srt, pt.size = .1, label = F, label.size = 4, group.by="compartment", reduction = "umap")
ggsave(paste0("compartments_UMAP.png"), width = 10, height = 8, type = "cairo")

message('Saving seurat object')
save(srt, file=paste0("srt.Rda"))


per_compartment <- SplitObject(srt, split.by = "compartment")

for (obj in per_compartment){
    name <- unique(obj$compartment)
    save(srt, file=paste0("srt_", name, ".Rda"))
}