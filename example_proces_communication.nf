process create_seurat_object {
    

    input:
        val(data_directory)

    output:
        path '*.csv', emit: out_csv

    script:
        """
        Rscript ${workflow.projectDir}/bin/1_create_seurat_object.R -d '$data_directory'
        """
}

dir = '/ahg/regevdata/projects/lungCancerBueno/Results/10x_nsclc_41421/data/PRIV_GITHUB/SCGA/nextflow/bin/data'


process second_test {
    

    input:
        path in_file

    output:
        stdout emit: echooo

    script:
        """
        Rscript ${workflow.projectDir}/bin/2_testing_file_save.R -csv '$in_file'
        """
}

workflow {    
    create_seurat_object(dir)
    second_test(create_seurat_object.out.flatten())
}


