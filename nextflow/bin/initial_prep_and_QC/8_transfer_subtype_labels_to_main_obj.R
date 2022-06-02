#!/usr/bin/env Rscript

library(argparse)
library(data.table)
library(ggplot2)
library(gridExtra)

parser <- ArgumentParser(description='An executible script to transfer subtype labels from compartment seurat objects back to the main seurat object.')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_object_path", help="Provide the full path to the seuart object .Rda file.")
parser$add_argument("-cso", "--compartment_seurat_objects", type="character", dest="cso_list", help="Provide a list of file paths with the final compartment seurat objects.")

args <- parser$parse_args()
seurat_object_path <- args$seurat_object_path
cso_list <- args$cso_list

message('Loading Main Seurat Object')
srt <- get(load(file=seurat_object_path))

# Load each of the individual seurat objects...
for (file_path in cso_list){
    assign(basename(file_path), get(load(file=file_path)))
}

# Get a vector of all seurat object names as strings
cso_names <- c()
for (file_path in cso_list){
    cso_names <- append(cso_names, basename(file_path))
}


# Use colnames and predicted.id column to map subtype labels back to the main object
srt$predicted.id <- 'placeholder'
for (name in cso_names){
    srt$predicted.id[colnames(srt) == colnames(eval(parse(text = name)))] <- eval(parse(text = paste0(name, '$predicted.id')))
}

# Need to make sure the placeholder is gone. It should be gone if the code ran correctly. I will be able to tell from a table which has the number of celltypes for each sampleID. It will be good to have the same table before and after doublet removal to see what cells have been removed.


meta_data <- seurat_object@meta.data
table <- table(meta_data[, c("sampleID", "compartment")])
table <- data.frame(table)
table <- dcast(table, sampleID ~ compartment)

#* Test this...
dcast(setDT(table), sampleID ~ compartment, value.var = c("predicted.id"))

# create output table
# Write to figure instead
png("cell_counts_per_sample_table.png")
p<-tableGrob(table)
grid.arrange(p)
dev.off()



