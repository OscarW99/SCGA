#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description='An executible script to re-cluster a compartment seurat object based on a given resolution.')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object_path", help="Provide the full path to the seuart object .Rda file.")
parser$add_argument("-pf", "--parameter_file", type="character", dest="parameter_file", help="Provide the path to a json file that contains clustering resolution")

#! The batch correction argument should be a string of all metadata columns seperated by a space. Also surround in quotes
args <- parser$parse_args()
seurat_object_path <- args$seurat_object_path
parameter_file <- args$parameter_file

message('Loading Seurat Object')
srt <- get(load(file=seurat_object_path))
# parameter file will have differnt name depending on compartment
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
#             "number_of_clusters":20,
#             "resolution":0.3
#         }
#     }
# }
# todo - check this
compartment <- unique(srt$compartment)[1]

message('Finding neighbours')
srt <- FindNeighbors(srt, dims = 1:30, verbose = FALSE)
message('Finding clusters')
resolution <- eval(parse(text = paste0("parameters$Compartments$",compartment,"$resolution")))
srt <- FindClusters(srt, verbose = FALSE, resolution = resolution)
Idents(srt) <- "seurat_clusters"

DimPlot(srt, pt.size = .1, label = F, label.size = 4, reduction = "umap")
ggsave(paste0(compartment, "_clusters_", resolution, "_UMAP.png"), width = 10, height = 8, type = "cairo")


message('Finding cluster markers')
#* 2) Top markers for each cluster in the form of a heatmap
srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ordered.srt.markers <- group_by(srt.markers, cluster)
top10 <- top_n(ordered.srt.markers, 10, wt = avg_log2FC)
DoHeatmap(srt, features = top10$gene, assay = "SCT") + NoLegend()
plot_size <- min(c(50, max(as.integer(srt$seurat_clusters))))
ggsave(paste0(compartment, "_cluster_markers_heatmap_", resolution, ".pdf"), width=(plot_size), height=(plot_size))
dev.off()

message('Plotting violin plots')
#* 4) A Violin-plot for each of the markers in the json file. This will show each clusters expresison of markers. Have all of these in one figure... so use cowplot.
input_markers <- unlist(eval(parse(text = paste0("parameters$Compartment_gene_sets$", compartment))))
vln_plots <- c()
for (i in 1:length(input_markers)){
    plot <- VlnPlot(srt, features=input_markers[[i]], log=TRUE, fill.by='feature',  pt.size=0) + 
    labs(title = paste0(input_markers[[i]] ," - ", names(input_markers)[i]))
    vln_plots[[i]] <- plot
}
if (length(vln_plots) < 8){
    fp_cols <- length(vln_plots)
    fp_width <- 3.5*length(vln_plots)
    fp_height <- 5
} else {
    fp_cols <- 8
    fp_width <- 28
    fp_height <- 3.5*((length(vln_plots) %% 8)+(length(vln_plots) %/% 8))
}
cowplot::plot_grid(plotlist = vln_plots, ncol = (fp_cols/2))
ggsave(paste0(compartment, "_violin_marker_plots", resolution, ".png"), width = fp_width, height = fp_height, type = "cairo")

parameters$Compartments[compartment]$number_of_clusters <- length(unique(srt.markers$cluster))

# Save the seurat object
message('Saving seurat object')
save(srt, file=paste0(compartment, "_srt_", resolution, ".Rda"))
write(parameters, paste0(compartment, "_", resolution, "_output.json"))



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
#             "resolution":0.3
#         }
#     }
# }