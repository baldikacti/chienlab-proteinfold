process AF3_FOLD {
    label "gpu_af3"
    publishDir "${params.outdir}/${params.mode}", mode: 'copy', pattern: "folds/*"

    input:
    path("input_a3m/*")
    path af3_db
    path af3_model

    output:
    path "folds/*"

    script:
    def args = task.ext.args ?: ''
    """
    module load alphafold3/latest
    mkdir folds
    run_alphafold.py --norun_data_pipeline --input_dir input_a3m --output_dir folds --db_dir $af3_db --model_dir $af3_model $args
    """
}