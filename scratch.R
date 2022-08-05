
library(ggplot2)
library(Seurat)
library(dplyr)

highlevel <- get(load(file="$PATH/highlevel.Rda"))



ggplot(highlevel@meta.data, aes(x=novelty_score)) + 
geom_density(alpha = 0.2) + 
#scale_x_log10() + 
theme_classic() +
ylab("Cell density") +
geom_vline(xintercept = 0.8, color='red')

ggsave("$PATH/test.png")



highlevel$novelty_score <- log10(length(rownames(highlevel)) / highlevel$nCount_RNA)