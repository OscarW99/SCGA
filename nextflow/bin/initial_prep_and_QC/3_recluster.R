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



message('Finding neighbours')
srt <- FindNeighbors(srt, dims = 1:30, verbose = FALSE)
message('Finding clusters')
resolution <- parameters$QC_thresholds$cluster_resolution
srt <- FindClusters(srt, verbose = FALSE, resolution = resolution)
Idents(srt) <- "seurat_clusters"

DimPlot(srt, pt.size = .1, label = F, label.size = 4, reduction = "umap")
ggsave(paste0("clusters_", resolution, "_UMAP.png"), width = 10, height = 8, type = "cairo")


message('Finding cluster markers')
#* 2) Top markers for each cluster in the form of a heatmap
srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ordered.srt.markers <- group_by(srt.markers, cluster)
top10 <- top_n(ordered.srt.markers, 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggsave(paste0("cluster_markers_heatmap", resolution, ".png"))

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
ggsave(paste0("Violin_marker_plots", resolution, ".png"), width = fp_width, height = fp_height, type = "cairo")

parameters$Clusters$number_of_clusters <- length(unique(srt$seurat_clusters))

# Save the seurat object
message('Saving seurat object')
save(srt, file=paste0("srt_", resolution, ".Rda"))
write(parameters, "output.json")