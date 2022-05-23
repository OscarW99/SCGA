#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description='An executible script to normalize, find variable features, scale, run PCA, run harmony batch correction and run (not plot) UMAPs')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object_path", help="Provide the full path to the seuart object .Rda file.")
parser$add_argument("-bc", "--batch_correct", type="character", dest="batch_correction", help="Provide the metadata columns of the seurat object that you want to perform batch correction on if any.")

#! The batch correction argument should be a string of all metadata columns seperated by a space. Also surround in quotes
args <- parser$parse_args()
seurat_object_path <- args$seurat_object_path
batch_correction <- args$batch_correction

if (is.null(batch_correction)){
    seurat_object <- NormalizeData(object = seurat_object)
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 3000)
    seurat_object <- ScaleData(seurat_object)
    seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
    seurat_object <- RunUMAP(seurat_object, reduction = "pca", dim = 1:15)
    seurat_object <- FindNeighbors(seurat_object, reduction = "pca", dims = 1:15)
    
} else {
    batch_correction <- unlist(strsplit(batch_correction, split=" "))
    batch_correction <- batch_correction[! batch_correction == ""]

    seurat_object <- NormalizeData(object = seurat_object)
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 3000)
    seurat_object <- ScaleData(seurat_object)
    seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
    seurat_object <- RunHarmony(seurat_object, batch_correction)
    seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dim = 1:15)
    seurat_object <- FindNeighbors(seurat_object, features = VariableFeatures(object = seurat_object))
}


save(seurat_object, file="normalized_UMAP_srt.Rda")

# todo - need some sort of error check, or rather spellcheck to see if the metadata columns provided for batch correction actually exist.