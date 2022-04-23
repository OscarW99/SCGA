#!/usr/bin/env nextflow


nextflow.enable.dsl=2




//@ CELLPHONEDB PROCESS SCRIPTS - WILL COMBINE INTO ONE WORKFLOW SCRIPT
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////

process cellphoneDB_split_seurat_object {
    
    
    input:
        // with only use path here if not using Rscript (val otherwise)
        tuple val(seurat_object_path), val(sample_id_meta)
        
    
    output:
        path '*rds', emit: patient_seurat_objects

    script:
        """
        Rscript ${workflow.projectDir}/bin/cellphoneDB/1_cellphoneDB_split_seurat_object.R -so $seurat_object_path -o ${workflow.worktDir} -id $sample_id_meta
        """
}

dir = "/ahg/regevdata/projects/lungCancerBueno/Results/10x_bischoff_102621/bischoff.obj.Rda"
sample_id_meta = "sampleID"

tuple = Channel.of([dir, sample_id_meta])
//file_prep_params = Channel.from[dir, sample_id_meta)


workflow {
    cellphoneDB_split_seurat_object(tuple)
    cellphoneDB_split_seurat_object.out.view()
}


process cellphoneDB_input_file_generation.R {
    
    
    input:
        // with only use path here if not using Rscript (val otherwise)
        tuple val(seurat_object_path), val(sample_id_meta)
        
    
    output:
        path '*rds', emit: patient_seurat_objects

    script:
        """
        Rscript ${workflow.projectDir}/bin/cellphoneDB/1_cellphoneDB_split_seurat_object.R -so $seurat_object_path -o ${workflow.worktDir} -id $sample_id_meta
        """
}


// process cellphoneDB_run {
    
//     // requires cpdb conda env
//     input:
//         path patient_count_data

//     output:
//         stdout emit: echooo

//     script:
//         """
//         Rscript ${workflow.projectDir}/bin/FIRST_TEST.R -so $seurat_object -bd $baseDir -c $compartment -n $draft
//         """
// }



// process cellphoneDB_dotplot {
    

//     input:
//         tuple val(compartment), val(draft), val(seurat_object), val(baseDir)

//     output:
//         stdout emit: echooo

//     script:
//         """
//         Rscript ${workflow.projectDir}/bin/FIRST_TEST.R -so $seurat_object -bd $baseDir -c $compartment -n $draft
//         """
// }




// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////

