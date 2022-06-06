#!/usr/bin/env nextflow


nextflow.enable.dsl=2



// process FIRST_TEST {
    

//     input:
//         tuple val(compartment), val(draft), val(seurat_object), val(baseDir)

//     output:
//         stdout emit: echooo

//     script:
//         """
//         Rscript ${workflow.projectDir}/bin/FIRST_TEST.R -so $seurat_object -bd $baseDir -c $compartment -n $draft
//         """
// }
// *Need to be careful not to change the names of scripts because once I name them it's a pain to go through the code as edit stuff.

process create_seurat_object {
    publishDir '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/PRIV_GITHUB/SCGA/nextflow/bin/publishDir', pattern: '*.png'

    input:
        val(data_directory)

    output:
        path 'srt.Rda', emit: seurat_out
        path 'nCount.png', emit: nCounts
        path 'nFeature.png', emit: nFeatures
        path 'percent.mt.png', emit: percent_mt


    script:
        """
        Rscript ${workflow.projectDir}/bin/1_create_seurat_object.R -d '$data_directory'
        """
}

dir = '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/PRIV_GITHUB/SCGA/nextflow/bin/data'


workflow {    
    create_seurat_object(dir)
    // second_test(create_seurat_object.out.flatten())
}





// process second_test {
    

//     input:
//         path in_file

//     output:
//         stdout emit: echooo

//     script:
//         """
//         Rscript ${workflow.projectDir}/bin/2_testing_file_save.R -csv '$in_file'
//         """
// }

// workflow {    
//     create_seurat_object(dir)
//     second_test(create_seurat_object.out.flatten())
// }










// FIRST_TEST
// highlevel_seurat_object = '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/pAdeno_early_naive_subset/htan_msk_addition/draft2/highlevel/v12.highlevel.Rda'
// baseDir = '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/pAdeno_early_naive_subset/htan_msk_addition/'
// draft = 'draft3'

// celltypes = Channel.from('bplasmast', 'myeloid', 'tnk')
// glob = Channel.of([draft, highlevel_seurat_object, baseDir])

// workflow {    
//     FIRST_TEST(celltypes.combine(glob))
//     FIRST_TEST.out.view()
//     celltypes.combine(glob).view()
// }