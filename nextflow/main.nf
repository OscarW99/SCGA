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
    publishDir '$PATH/SCGA/nextflow/bin/publishDir', pattern: '*.png'

    input:
        val(data_directory)

    output:
        path 'srt.Rda', emit: seurat_out
        path 'nCount.png', emit: nCounts
        path 'nFeature.png', emit: nFeatures
        path 'percent.mt.png', emit: percent_mt


    script:
        """
        Rscript ${workflow.projectDir}/bin/initial_prep_and_QC/1_create_seurat_object.R -d '$data_directory'
        """
}

dir = '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/PRIV_GITHUB/SCGA/nextflow/bin/data'

// // todo - test this script
process qc_filtering {
    publishDir '$PATH/SCGA/nextflow/bin/publishDir', pattern: '*.{png/pdf}'

    input:
        tuple val(seurat_object_path), val(parameter_file)

    output:
        path '*.png', emit: png_figs
        path '*.pdf', emit: pdf_figs
        path 'output.json', emit: json_out
        path '*.Rda', emit: seurat_out

    script:
        """
        Rscript ${workflow.projectDir}/bin/initial_prep_and_QC/2_qc_filtering.R -so '$seurat_object_path' -pf '$parameter_file'
        """
}

dir2 = '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/PRIV_GITHUB/SCGA/2_inputs.json'

srt = '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/PRIV_GITHUB/SCGA/nextflow/work/1c/dab97150c674743f2156737e04fa19/srt.Rda'

workflow {   
    // create_seurat_object(dir)
    // qc_in = Channel.of( [create_seurat_object.out.seurat_out, dir2] )
    qc_in = Channel.of( [srt, dir2] )
    qc_filtering( qc_in )
}




// FIRST_TEST
// highlevel_seurat_object = '$PATH/v12.highlevel.Rda'
// baseDir = '$PATH/base/'
// draft = 'draft3'

// celltypes = Channel.from('bplasmast', 'myeloid', 'tnk')
// glob = Channel.of([draft, highlevel_seurat_object, baseDir])

// workflow {    
//     FIRST_TEST(celltypes.combine(glob))
//     FIRST_TEST.out.view()
//     celltypes.combine(glob).view()
// }