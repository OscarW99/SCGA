
library(ggplot2)
library(Seurat)
library(dplyr)

highlevel <- get(load(file="/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/pAdeno_early_naive_subset/htan_msk_addition/draft2/highlevel/with_subtypes.highlevel.Rda"))



ggplot(highlevel@meta.data, aes(x=novelty_score)) + 
geom_density(alpha = 0.2) + 
#scale_x_log10() + 
theme_classic() +
ylab("Cell density") +
geom_vline(xintercept = 0.8, color='red')

ggsave("/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/PRIV_GITHUB/figs/test.png")



highlevel$novelty_score <- log10(length(rownames(highlevel)) / highlevel$nCount_RNA)