#!/usr/bin/env Rscript


library(argparse)

parser <- ArgumentParser(description='An executible script to prepare files to run cellphoneDB in a patient-wise manner.')

parser$add_argument("-d", "--data_directory", type="character", dest="data_directory", help="Provide the full path to directory in which all patient count data is stored.")

args <- parser$parse_args()
data_directory <- args$data_directory

#! del these
data_directory <- '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/PRIV_GITHUB/SCGA/nextflow/bin/data'
#!


















library(ggplot2)
library(Seurat)
library(dplyr)
library(Matrix)

data.path <- "/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/pAdeno_early_naive_subset/htan_msk_addition/draft2/highlevel/v12.highlevel.Rda"
proj.path <- "/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/cellphonedb/new_age"


# source(paste0(user.path, "/R_utils/plotutils.R"))
# source(paste0(user.path, "/R_utils/seuratutils.R"))
# source(paste0(user.path, "/R_utils/seuratutilsV3.R"))
# source(paste0(user.path, "/R_utils/color.R"))

# user.path <- "/ahg/regevdata/projects/ICA_Lung/Will"
# #Genesets
# cc.genes <- readLines(con = paste0(user.path, "/genelists/regev_lab_cell_cycle_genes.txt"))
# s.genes <- cc.genes[1:43]
# g2m.genes <- cc.genes[44:97]

# figures.dir <- paste0(proj.path, "/Results/10x_All_190615/Allv9/")
# load(file = paste0(figures.dir, "all.merge.new.cci.nodoublets.v2.Rda"))
load(file = data.path)
# all.bypatient <- SplitObject(all, split.by = "orig.identSec")
# split by patientID
all.bypatient <- SplitObject(v12.celltype, split.by = "sampleID")

ls()
typeof(all.bypatient)

library(future.apply)

cpdb.output.path <- paste0(proj.path, "/cell_cell_int_padeno_highlvl/")
dir.create(file.path(cpdb.output.path), showWarnings = FALSE)

options(future.globals.maxSize = 48000 * 1024^2)



#*############
saveRDS(all.bypatient, file = paste0(proj.path, "/all.bypatient.rds"))

sample.names <- list.files("/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/cellphonedb/")
sample.names <- sample.names[!(sample.names %in% c("draft_results","smoking_high_low.csv", "Differential_analysis", "new_Differential_analysis","EGFR_mut_wt.csv", "new_age"))]
old_samples <- vector("numeric")
for (i in 1:length(all.bypatient)){
  if (unique(all.bypatient[[i]]$sampleID) %in% sample.names){
    print(unique(all.bypatient[[i]]$sampleID))
    old_samples <- append(old_samples, i)
  }
}
print(sample.names)

# Only want to run cellphonedb on samples we haven't run before.
new.bypatient <- all.bypatient[-old_samples]
saveRDS(new.bypatient, file = paste0(proj.path, "/new.bypatient.rds"))


#*############


data.path <- "/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/pAdeno_early_naive_subset/htan_msk_addition/draft2/highlevel/v12.highlevel.Rda"
proj.path <- "/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/cellphonedb/new_age"

new.bypatient <- readRDS(paste0(proj.path, "/new.bypatient.rds"))

library(future.apply)
cpdb.output.path <- paste0(proj.path, "/cell_cell_int_padeno_highlvl/")
options(future.globals.maxSize = 48000 * 1024^2)

sub_new.bypatient <- new.bypatient[31:35] #! ... EDIT THIS


lapply(sub_new.bypatient, function(x){
  sampleid <- unique(x$sampleID)
  print(sampleid)
  print(paste0("Creating directory:", paste0(cpdb.output.path, sampleid)))
  dir.create(file.path(cpdb.output.path, sampleid), showWarnings = FALSE)
  counts <- x[["RNA"]]@counts
  counts.norm <- future_apply(counts, 2, function(x) (x/sum(x))*10000) # this is recommended by the cellphonedb paper
  print('2')
  write.table(counts.norm, paste0(cpdb.output.path, sampleid, "/", sampleid, "_counts_highlvl.txt"), sep="\t", quote=F)
  print(paste0(cpdb.output.path, sampleid, "/", sampleid, "_counts_highlvl.txt"))
  meta <- x@meta.data
  meta <- cbind(rownames(meta), meta[,'celltype_highlevel', drop=F])
  meta[is.na(meta[,2]),2]<-'Unassigned'
  colnames(meta)<-c('cell','annotation')
  write.table(meta, paste0(cpdb.output.path, sampleid, "/", sampleid, "_meta_highlvl.txt"), sep="\t", quote=F, row.names=F)
  print(paste0(cpdb.output.path, sampleid, "/", sampleid, "_meta_highlvl.txt"))
})