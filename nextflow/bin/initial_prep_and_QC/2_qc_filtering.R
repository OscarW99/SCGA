#!/usr/bin/env Rscript#!/usr/bin/env Rscript


library(argparse)
library(Seurat)
library(rjson)
library(ggplot2)
library(sctransform)
library(dplyr)

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
#         "Myeloid":["S100A8", "RETN", "FCN1"],
#         "Bplasmast":["IGHGP", "IGHG1", "IGHG3"]
#     }
# }

cells_before_filter <- length(colnames(srt))

# QC filtering
srt <- subset(srt, subset = nFeature_RNA > as.integer(parameters$QC_thresholds$nFeatures_threshold) & nCount_RNA > as.integer(parameters$QC_thresholds$nCounts_threshold) & percent.mt < as.integer(parameters$QC_thresholds$percent_mt_threshold))

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
srt <- FindClusters(srt, verbose = FALSE, resolution = as.integer(parameters$QC_thresholds$cluster_resolution)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
srt <- CellCycleScoring(srt, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)


#* 1) UMAPs for clusters, cell cycle score, %mt, nCount, nFeatures
DimPlot(srt, pt.size = .1, label = F, label.size = 4, reduction = "umap")
ggsave("clusters_UMAP.png"), width = 10, height = 8, type = "cairo")
# !##############################################################################
DimPlot(srt, pt.size = .1, label = F, label.size = 4, group.by = "Phase", reduction = "umap")
ggsave("cellcycle_phase_UMAP.png"), width = 10, height = 8, type = "cairo")
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


#* 2) Top markers for each cluster in the form of a heatmap
srt.markers <- FindAllMarkers(srt, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ordered.srt.markers <- group_by(srt.markers, cluster)
top10 <- top_n(ordered.srt.markers, 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
ggsave("cluster_markers_heatmap.png")


#* 3) A Feature-plot with all celltype markers. I could paste the celltype before the gene name for the title of each feature plot.
# Get features by looping through all json markers
input_markers <- unlist(parameters$Gene_sets)

p <- FeaturePlot(seurat_object, input_markers, pt.size = .001, combine = FALSE)
for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
lapply(p, function(x){x + labs(title = names(input_markers))})
cowplot::plot_grid(plotlist = p, ncol = 8)
ggsave("Feature_plot.png", type = "cairo")

# TODO - test the above works







#* 4) A Violin-plot for each of the markers in the json file. This will show each clusters expresison. Maybe produce pdf with all violin plots... or can I somehow show multiple plots in Rshiny in a optimised way?

# 3) I will want to be able to play around with different plots etc. I think saving the srt object here and then 
#  The next script will be where I assign celltpye names to the cells


#* COMMENTS -------------------------------------------------------------------------------------------------
# I should do somthing like those people at mt-sinai whos nextflow pipeline produces pdf slides with each figure explaining what they are. I could also then include what the statistics show ect.

#! How will I Output the number of cells filtered? As a graph, A table? both. Also is there an R package to create a publication quality figure. It would be good to auto-generate figures and figure legnends.

#? Based on the tissue type I could have known markers plotted with markers in this tsv file... https://panglaodb.se/markers.html?cell_type=%27choose%27
#? I may have to re-download this file periodically for updated versions? Although it hasn't been updated for 2 years so perhaps I wont need to do this.
#? I can also add my own markers to this file. So at least for the lung I can add them and then perhaps manually adding ones per tissue based on the literatire for specific clients. Eventually building up a marker database perhaps for mouse and human.

# Todo - print outputs to the console so I can see whats going on, i.e what stage in this script is currently being run.

# Need seperate nf script to recluster srt at a different resolution... and then re-run all plots.