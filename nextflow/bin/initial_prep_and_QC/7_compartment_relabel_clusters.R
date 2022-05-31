#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description='An executible script to re-cluster a compartment seurat object based on a given resolution.')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object_path", help="Provide the full path to the seuart object .Rda file.")
parser$add_argument("-pf", "--parameter_file", type="character", dest="parameter_file", help="Provide the path to a json file that contains clustering resolution")

args <- parser$parse_args()
seurat_object_path <- args$seurat_object_path
parameter_file <- args$parameter_file

message('Loading Seurat Object')
srt <- get(load(file=seurat_object_path))
parameters <- fromJSON(file=parameter_file)

# What json will look like at the end
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
#     "Compartment_gene_sets":{
#         "Myeloid":{
#             "Momac":["CD80","HLA-DR", "CD206"],
#             "CD14.Mono":["xyz", "abc", "lmn"]
#         },
#         "Bplasmast":{
#             "first_bst":["b_one","b_two", "b_three"],
#             "second_bst":["b_a", "b_b", "b_c"]
#         }
#     },
#     "Clusters":{
#         "number_of_clusters":5
#         "new_cluster_names":c("Tcell", "Bcel", "Myeloid", "Epithelial", "Endothelial")
#     },
#     "Compartments":{
#         "Myeloid":{
#             "number_of_clusters":16,
#             "resolution":0.3,
#             "new_cluster_names":c("moMac_APOE", "DC1", "DC2", "cycling_cell", "monocytes", "...")
#         }
#     }
# }
compartment <- unique(srt$compartment)[1]

new_labels <- parameters$Compartments[compartment]$new_cluster_names
srt$predicted.id <- lapply(srt$seurat_clusters, function(x) new_labels[x+1])

DimPlot(srt, pt.size = .1, label = F, label.size = 4, group.by="predicted.id", reduction = "umap")
ggsave(paste0(compartment, "_compartments_UMAP.png"), width = 10, height = 8, type = "cairo")

message('Saving seurat object')
save(srt, file=paste0(compartment, ".Rda"))
