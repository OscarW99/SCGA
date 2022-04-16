#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process low_level_integration {
    

    input:
        tuple val(compartment), val(draft), val(seurat_object), val(baseDir)

    output:
        stdout emit: echooo

    script:
        """
        Rscript ${workflow.projectDir}/bin/low_level_label_transfer.R -so $seurat_object -bd $baseDir -c $compartment -n $draft
        """

}


highlevel_seurat_object = '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/pAdeno_early_naive_subset/htan_msk_addition/draft2/highlevel/v12.highlevel.Rda'
baseDir = '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/pAdeno_early_naive_subset/htan_msk_addition/'
draft = 'draft3'

celltypes = Channel.from('bplasmast', 'myeloid', 'tnk')
glob = Channel.of([draft, highlevel_seurat_object, baseDir])

workflow {    
    low_level_integration(celltypes.combine(glob))
    low_level_integration.out.view()
}