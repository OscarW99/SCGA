#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser(description='An executible script to run compartment-wise label transfer. Takes x arguments, a compartment, ')

parser$add_argument("-so", "--seurat_object", type="character", dest="seurat_path", help="Provide the full path to the highlevel seurat object.")

parser$add_argument("-bd", "--baseDir", type="character", dest="baseDir", help="Provide the full path to the directory that you want the output to be stored.")

parser$add_argument("-c", "--compartment", type="character", dest="compartment", help="Provid the name of the compartment in your seurat object.")

parser$add_argument("-n", "--name", type="character", dest="name", help="Provid a name for this iteration of compartment-wise label transfer. E.g. 'draft1'")

args <- parser$parse_args()
print(args)

source("/home/unix/owright/utils/base_functions.R")

print(getwd())
print('Done')