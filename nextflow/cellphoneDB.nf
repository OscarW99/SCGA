#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//@ CELLPHONEDB PROCESS SCRIPTS - WILL COMBINE INTO ONE WORKFLOW SCRIPT
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////

// //* Which one is it?
// projectDir
// launchDir
// workDir


process cellphoneDB_split_seurat_object {
    
    
    input:
        // with only use path here if not using Rscript (val otherwise)
        tuple val(seurat_object_path), val(sample_id_meta), val(celltype_label_meta)
        
    
    output:
        //path '*rds', emit: patient_seurat_objects
        // might need to use a star * instead of {x} below
        path "*_meta.txt", emit: meta_file
        path "*_counts.txt", emit: counts_file
        // stdout emit: verbiage

    script:
        """
        Rscript ${workflow.projectDir}/bin/cellphoneDB_scripts/1_cellphoneDB_split_seurat_object.R -so $seurat_object_path -id $sample_id_meta -l $celltype_label_meta
        """
}

dir = "/ahg/regevdata/projects/lungCancerBueno/Results/10x_bischoff_102621/bischoff.obj.Rda"
sample_id_meta = "sampleID"
celltype_label_meta = "predicted.id.highlevel"

tuple = Channel.of([dir, sample_id_meta, celltype_label_meta])



// cpdb_file_prep
//     .fromFilePairs('${workDir}/*_{meta,counts}_highlvl.txt')


process cellphoneDB_run {
    conda "/ahg/regevdata/projects/ICA_Lung/Oscar/conda/cpdb"
    // cpus 1
    // executor 'sge'
    // memory '24 GB'
    // time '5h'


    input:
        tuple val(counts_file), val(meta_file)

    output:
        stdout emit: echooo
        // cellphonedb method statistical_analysis --counts-data=gene_name $meta_file	$counts_file

    script:
        """
        cat $counts_file $meta_file > file4.txt
        """
}

workflow {
    cellphoneDB_split_seurat_object(tuple)
    cellphoneDB_run(cellphoneDB_split_seurat_object.out.collect())
    cellphoneDB_run.out.view()
}


// cpdb_files = Channel.fromPath( '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/adj_normal_subset/for_manu/draft3/highlevel/cell_cell_int_prep/*')
// cpdb_files.flatten().view()

// workflow {
//     cellphoneDB_run(cpdb_files.flatten())
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































// // -----------------------------------------

// process cellphoneDB_input_file_generation {
    
    
//     input:
//         // with only use path here if not using Rscript (val otherwise)
//         tuple val(seurat_object_path), val(sample_id_meta), val(celltype_label_meta)
        
    
//     output:
//         stdout emit: echooo

//     script:
//         """
//         Rscript ${workflow.projectDir}/bin/cellphoneDB_scripts/2_cellphoneDB_input_file_generation.R -so $seurat_object_path -o ${workflow.workDir} -id $sample_id_meta -l $celltype_label_meta
//         """
// }


// seurat_objects = Channel.fromPath( '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/adj_normal_subset/for_manu/draft3/highlevel/cellphoneDB/patient_seurat_objects/*.rds')


// sample_id_meta = "SampleID"
// celltype_label_meta = "luad_label_match"
// constants = Channel.of([sample_id_meta, celltype_label_meta])

// // // // //* I need to create a channel from the dir of the part1 output.


// workflow {
//     // cellphoneDB_split_seurat_object(tuple)
//     // seurat_objects = Channel.from(cellphoneDB_split_seurat_object.out.flatten())
//     cellphoneDB_input_file_generation(seurat_objects.combine(constants))
//     cellphoneDB_input_file_generation.out.view()
    
    
// }

// // -----------------------------------------









// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////
// @ //////////////////////////////////////////////////////////////////////////////////////////

