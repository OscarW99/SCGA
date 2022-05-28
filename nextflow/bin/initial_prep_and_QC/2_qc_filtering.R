#!/usr/bin/env Rscript#!/usr/bin/env Rscript

# todo - can i somehow automate testing, like with pytest. I may make a lot of chamges to these scripts going forward so it would be a good idea to make tests. ATM I am just loading stuff in from random places so I could atleast have commented out tests.
library(argparse)
library(Seurat)
library(rjson)
library(ggplot2)
library(sctransform)
library(dplyr)
library(gridExtra)

parser <- ArgumentParser(description='An executible R script to filter cells in a seurta object based on provided thresholds for nCounts, nFeatures and percent.mt genes.')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object_path", help="Provide the full path to the seuart object .Rda file.")
parser$add_argument("-pf", "--parameter_file", type="character", dest="parameter_file", help="Provide the path to a json file that contains threshold values for nCounts, nFeatures and percent.mt genes.")

args <- parser$parse_args()
seurat_object_path <- args$seurat_object_path
parameter_file <- args$parameter_file

message('Loading Seurat Object')

# TESTING
# seurat_object_path <- '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/pAdeno_early_naive_subset/htan_msk_addition/draft2/myeloid/v12.myeloid.Rda'
# parameter_file <- '2_inputs.json'
###
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
#         "Myeloid":["S100A8", "RETN", "FCN1"],
#         "Bplasmast":["IGHGP", "IGHG1", "IGHG3"]
#     }
# }

cells_before_filter <- length(colnames(srt))

# QC filtering
message('Subsetting Seurat Object')
srt <- subset(srt, subset = nFeature_RNA > as.integer(parameters$QC_thresholds$nFeatures_threshold) & nCount_RNA > as.integer(parameters$QC_thresholds$nCounts_threshold) & percent.mt < as.integer(parameters$QC_thresholds$percent_mt_threshold))

cells_after_filter <- length(colnames(srt))

# Regressing out the difference between the G2M and S phase scores. This means that signals separating non-cycling cells and cycling cells will be maintained, but differences in cell cycle phase amongst proliferating cells (which are often uninteresting), will be regressed out of the data
# marrow$CC.Difference <- marrow$S.Score - marrow$G2M.Score
# marrow <- ScaleData(marrow, vars.to.regress = "CC.Difference", features = rownames(marrow))
# marrow <- RunPCA(marrow, features = VariableFeatures(marrow), nfeatures.print = 10)


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
resolution <- parameters$QC_thresholds$cluster_resolution
srt <- FindClusters(srt, verbose = FALSE, resolution = resolution)
Idents(srt) <- "seurat_clusters"

message('Performing cell cycle analysis')
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
srt <- CellCycleScoring(srt, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

message('Plotting UMAPs')
#* 1) UMAPs for clusters, cell cycle score, %mt, nCount, nFeatures
DimPlot(srt, pt.size = .1, label = F, label.size = 4, reduction = "umap")
ggsave(paste0("clusters_", resolution, "UMAP.png"), width = 10, height = 8, type = "cairo")
# !##############################################################################
DimPlot(srt, pt.size = .1, label = F, label.size = 4, group.by = "Phase", reduction = "umap")
ggsave("cellcycle_phase_UMAP.png", width = 10, height = 8, type = "cairo")
# !##############################################################################
n_Feature_df <- data.frame(srt@reductions$umap@cell.embeddings, nFeature_RNA = srt@meta.data$nFeature_RNA)
p <- ggplot(data = n_Feature_df) +
    geom_point(mapping = aes(UMAP_1, UMAP_2, color = log10(nFeature_RNA)), size = 0.01) +
    scale_colour_gradientn(colours = rainbow(7))
ggsave("nFeature_UMAP.png", width = 10, height = 5, type = "cairo")
# !##############################################################################
n_count_df <- data.frame(srt@reductions$umap@cell.embeddings, nCount_RNA = srt@meta.data$nCount_RNA)
p <- ggplot(data = n_count_df) +
    geom_point(mapping = aes(UMAP_1, UMAP_2, color = log10(nCount_RNA)), size = 0.01) +
    scale_colour_gradientn(colours = rainbow(7), limits = quantile(log10(n_count_df$nCount_RNA), c(0, 1)))
ggsave("nCount_UMAP.png", width = 10, height = 5, type = "cairo")
# !##############################################################################
percmito_df <- data.frame(srt@reductions$umap@cell.embeddings, percent_mitochondrial_genes = srt@meta.data$percent.mt)
p <- ggplot(data = percmito_df) +
    geom_point(mapping = aes(UMAP_1, UMAP_2, color = percent_mitochondrial_genes), size = 0.01) +
    scale_colour_gradientn(colours = rainbow(7))
ggsave("percent_mitochondrial_UMAP.png", width = 10, height = 5, type = "cairo")
# !##############################################################################

message('Finding cluster markers')
#* 2) Top markers for each cluster in the form of a heatmap
srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ordered.srt.markers <- group_by(srt.markers, cluster)
top10 <- top_n(ordered.srt.markers, 10, wt = avg_log2FC)
DoHeatmap(srt, features = top10$gene) + NoLegend()
ggsave(paste0("cluster_markers_heatmap", resolution, ".png"))

message('Plotting feature plots')
#* 3) A Feature-plot with all celltype markers. I could paste the celltype before the gene name for the title of each feature plot.
input_markers <- unlist(parameters$Gene_sets)

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
ggsave("Feature_plot.png", width = fp_width, height = fp_height, type = "cairo")


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


# Create table of filtered cells and output this.
Stage <- c("Before filtering", "After filtering", "number removed")
Cells <- c(cells_before_filter, cells_after_filter, (cells_before_filter-cells_after_filter))
filtering_df <- data.frame(Stage, Cells)
png("cells_filtered.png")
t<-tableGrob(filtering_df)
grid.arrange(t)
dev.off()

# update parameters json with output info
parameters$Cell_filtering$cells_before_filtering <- cells_before_filter
parameters$Cell_filtering$cells_after_filtering <- cells_after_filter
parameters$Cell_filtering$cells_filtered <- (cells_before_filter-cells_after_filter)

parameters$Clusters$number_of_clusters <- length(unique(srt.markers$cluster))


# Save the seurat object
message('Saving seurat object')
save(srt, file=paste0("srt_", resolution, ".Rda"))
write(parameters, "output.json")

# *# the json output will be the same for any clustering resolution. Only the seuratobject file name/ contents will change.
#* COMMENTS -------------------------------------------------------------------------------------------------
# I should do somthing like those people at mt-sinai whos nextflow pipeline produces pdf slides with each figure explaining what they are. I could also then include what the statistics show ect.

#! How will I Output the number of cells filtered? As a graph, A table? both. Also is there an R package to create a publication quality figure. It would be good to auto-generate figures and figure legnends.

#? Based on the tissue type I could have known markers plotted with markers in this tsv file... https://panglaodb.se/markers.html?cell_type=%27choose%27
#? I may have to re-download this file periodically for updated versions? Although it hasn't been updated for 2 years so perhaps I wont need to do this.
#? I can also add my own markers to this file. So at least for the lung I can add them and then perhaps manually adding ones per tissue based on the literatire for specific clients. Eventually building up a marker database perhaps for mouse and human.

# Need seperate nf script to recluster srt at a different resolution... and then re-run all plots.