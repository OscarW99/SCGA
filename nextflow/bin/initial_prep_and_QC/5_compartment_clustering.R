
library(argparse)
library(Seurat)
library(rjson)
library(ggplot2)
library(sctransform)
library(dplyr)
library(gridExtra)

parser <- ArgumentParser(description='An executible R script to perform de novo clustering per compartment.')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object_path", help="Provide the full path to the seuart object .Rda file.")
parser$add_argument("-pf", "--parameter_file", type="character", dest="parameter_file", help="Provide the path to a json file that contains compartment subtype markers.")

args <- parser$parse_args()
seurat_object_path <- args$seurat_object_path
parameter_file <- args$parameter_file

message('Loading Seurat Object')
srt <- get(load(file=seurat_object_path))
parameters <- fromJSON(file=parameter_file)


# # What json import will look like
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
#     }
# }


options(future.globals.maxSize = 4000 * 1024^2)
message('Performing SCTransform')
srt <- SCTransform(srt, vars.to.regress = "percent.mt", verbose = FALSE)
message('Running PCA')
srt <- RunPCA(srt, verbose = FALSE)
message('Running UMAP')
srt <- RunUMAP(srt, dims = 1:30, verbose = FALSE)
message('Finding neighbours')
srt <- FindNeighbors(srt, dims = 1:30, verbose = FALSE)
message('Finding clusters')
resolution <- 0.5
srt <- FindClusters(srt, verbose = FALSE, resolution = resolution)
Idents(srt) <- "seurat_clusters"

compartment <- unique(srt$compartment)[1]


##*#################################################

message('Plotting UMAPs')
#* 1) UMAPs for clusters, cell cycle score, %mt, nCount, nFeatures
DimPlot(srt, pt.size = .1, label = F, label.size = 4, reduction = "umap")
ggsave(paste0(compartment, "_clusters_", resolution, "UMAP.png"), width = 10, height = 8, type = "cairo")
# !##############################################################################
DimPlot(srt, pt.size = .1, label = F, label.size = 4, group.by = "Phase", reduction = "umap")
ggsave(paste0(compartment, "_cellcycle_phase_UMAP.png"), width = 10, height = 8, type = "cairo")
# !##############################################################################
n_Feature_df <- data.frame(srt@reductions$umap@cell.embeddings, nFeature_RNA = srt@meta.data$nFeature_RNA)
p <- ggplot(data = n_Feature_df) +
    geom_point(mapping = aes(UMAP_1, UMAP_2, color = log10(nFeature_RNA)), size = 0.01) +
    scale_colour_gradientn(colours = rainbow(7))
ggsave(paste0(compartment, "_nFeature_UMAP.png"), width = 10, height = 5, type = "cairo")
# !##############################################################################
n_count_df <- data.frame(srt@reductions$umap@cell.embeddings, nCount_RNA = srt@meta.data$nCount_RNA)
p <- ggplot(data = n_count_df) +
    geom_point(mapping = aes(UMAP_1, UMAP_2, color = log10(nCount_RNA)), size = 0.01) +
    scale_colour_gradientn(colours = rainbow(7), limits = quantile(log10(n_count_df$nCount_RNA), c(0, 1)))
ggsave(paste0(compartment, "_nCount_UMAP.png"), width = 10, height = 5, type = "cairo")
# !##############################################################################
percmito_df <- data.frame(srt@reductions$umap@cell.embeddings, percent_mitochondrial_genes = srt@meta.data$percent.mt)
p <- ggplot(data = percmito_df) +
    geom_point(mapping = aes(UMAP_1, UMAP_2, color = percent_mitochondrial_genes), size = 0.01) +
    scale_colour_gradientn(colours = rainbow(7))
ggsave(paste0(compartment, "_percent_mitochondrial_UMAP.png"), width = 10, height = 5, type = "cairo")
# !##############################################################################

message('Finding cluster markers')
#* 2) Top markers for each cluster in the form of a heatmap
srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ordered.srt.markers <- group_by(srt.markers, cluster)
top10 <- top_n(ordered.srt.markers, 10, wt = avg_log2FC)
DoHeatmap(srt, features = top10$gene) + NoLegend()
ggsave(paste0(compartment, "_cluster_markers_heatmap", resolution, ".png"))

message('Plotting feature plots')
#* 3) A Feature-plot with all celltype markers. I could paste the celltype before the gene name for the title of each feature plot.

# compartment is case sensitive
# todo - test this
input_markers <- unlist(eval(parse(text = paste0("parameters$Gene_sets$", compartment))))

p <- FeaturePlot(srt, input_markers, pt.size = .001, combine = FALSE)
for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + NoLegend() + NoAxes() + labs(title = paste0(input_markers[[i]] ," - ", names(input_markers)[i]))
}
if (length(p) < 8){
    fp_cols <- length(p)
    fp_width <- 3.5*length(p)
    fp_height <- 5
} else {
    fp_cols <- 8
    fp_width <- 28
    fp_height <- 3.5*((length(p) %% 8)+(length(p) %/% 8))
}
cowplot::plot_grid(plotlist = p, ncol = fp_cols)
ggsave(paste0(compartment, "_feature_plot.png"), width = fp_width, height = fp_height, type = "cairo")


message('Plotting violin plots')
#* 4) A Violin-plot for each of the markers in the json file. This will show each clusters expresison of markers. Have all of these in one figure... so use cowplot.
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

parameters$Compartments[compartment] <- length(unique(srt.markers$cluster))


# Save the seurat object
message('Saving seurat object')
save(srt, file=paste0(compartment, "_srt_", resolution, ".Rda"))
write(parameters, paste0(compartment, "_output.json"))# rename after compartment

# # What json import will look like
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
#     }
# }










# script 6 (reclustering needs to have new compartment clustering resolutions. These must be retrieved from the Rshiny front end.)