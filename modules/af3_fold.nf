process AF3_FOLD {
    label 'gpu'
    label 'error_ignore'
    publishDir "${params.outdir}/${params.mode}", mode: 'copy', pattern: "folds/*"

    container "docker://baldikacti/alphafold3:latest"

    input:
    path("input_a3m/*")
    path af3_db
    path af3_model

    output:
    path "folds/*"
    path ("folds/*/*_summary_confidences.json") , emit: summary_json

    script:
    def args = task.ext.args ?: ''
    """
    mkdir folds
    run_alphafold.py --norun_data_pipeline --input_dir input_a3m --output_dir folds --db_dir $af3_db --model_dir $af3_model $args
    """
}