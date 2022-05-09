#!/usr/bin/env Rscript


library(argparse)

parser <- ArgumentParser(description='Testing file finder')

parser$add_argument("-csv", "csv-file", type="character", dest="file", help="Path to csv file")

args <- parser$parse_args()
file <- args$file



in_file <- read.csv(file = file)

write.csv(in_file)
write.csv("hello mate")