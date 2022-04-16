#!/usr/bin/env Rscript

# For this one we will start by assuming the data will be in the form of count data for each patient and the folders will be named after the sampleID.

# So using Bischoff as an example for now. We wiil need to see what the actual output from 10X ect is. We could design the pipeline to integrate seamlessly from these sequencing technologies.


## EXAMPLE DIRECOTRY STRUCTURE ##
# |-- patient1
# |   `-- filtered_feature_bc_matrix
# |       |-- barcodes.tsv.gz
# |       |-- features.tsv.gz
# |       `-- matrix.mtx.gz
# |-- patient2
# |   `-- filtered_feature_bc_matrix
# |       |-- barcodes.tsv.gz
# |       |-- features.tsv.gz
# |       `-- matrix.mtx.gz
# `-- patient2
#     `-- filtered_feature_bc_matrix
#         |-- barcodes.tsv.gz
#         |-- features.tsv.gz
#         `-- matrix.mtx.gz


## Arguments?
# The root folder for these files
library(argparse)

parser <- ArgumentParser(description='An executible script to create a merged seurat object from the single cell sequencing data from many patients.')

parser$add_argument("-d", "--data_directory", type="character", dest="data_directory", help="Provide the full path to directory in which all patient count data is stored.")

args <- parser$parse_args()
data_directory <- args$data_directory

# print(args)


# print(getwd())
patient_folders <- list.dirs(data_directory)

print(patient_folders)
