#!/usr/bin/env nextflow


nextflow.enable.dsl=2




//@ CELLPHONEDB PROCESS SCRIPTS - WILL COMBINE INTO ONE WORKFLOW SCRIPT
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////

process cellphoneDB_split_seurat_object {
    
    
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

