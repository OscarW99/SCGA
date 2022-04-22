#* rename multiple seurat columns at once...
metadata <- merged_seurat@meta.data
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                    nUMI = nCount_RNA,
                    nGene = nFeature_RNA)