#!/usr/bin/env nextflow


nextflow.enable.dsl=2


//@ CELLPHONEDB PROCESS SCRIPTS - WILL COMBINE INTO ONE WORKFLOW SCRIPT
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////

process cellphoneDB_input_file_generation {
    

    input:
        // with only use path here if not using Rscript (val otherwise)
        val seurat_object_path
        

    output:
        path or val?

    script:
        """
        Rscript ${workflow.projectDir}/bin/cellphoneDB_input_file_generation.R -so $seurat_object_path -o ${workflow.worktDir}
        """
}




process cellphoneDB_run {
    
    // requires cpdb conda env
    input:
        path patient_count_data

    output:
        stdout emit: echooo

    script:
        """
        Rscript ${workflow.projectDir}/bin/FIRST_TEST.R -so $seurat_object -bd $baseDir -c $compartment -n $draft
        """
}



process cellphoneDB_dotplot {
    

    input:
        tuple val(compartment), val(draft), val(seurat_object), val(baseDir)

    output:
        stdout emit: echooo

    script:
        """
        Rscript ${workflow.projectDir}/bin/FIRST_TEST.R -so $seurat_object -bd $baseDir -c $compartment -n $draft
        """
}




// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////



process FIRST_TEST {
    

    input:
        tuple val(compartment), val(draft), val(seurat_object), val(baseDir)

    output:
        stdout emit: echooo

    script:
        """
        Rscript ${workflow.projectDir}/bin/FIRST_TEST.R -so $seurat_object -bd $baseDir -c $compartment -n $draft
        """
}
// *Need to be careful not to change the names of scripts because once I name them it's a pain to go through the code as edit stuff.

process create_seurat_object {
    

    input:
        val(data_directory)

    output:
        stdout emit: echooo

    script:
        """
        Rscript ${workflow.projectDir}/bin/1_create_seurat_object.R -d '$data_directory'
        """
}

dir = "/ahg/regevdata/projects/lungCancerBueno/Results/10x_bischoff_102621/data/"

workflow {    
    create_seurat_object(dir)
    create_seurat_object.out.view()
}










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