process AF3_FOLD {
    label 'process_medium'
    // toString() guards against a command-line `--use_gpu false`, which reaches a
    // process directive as the String 'false' rather than a boolean
    label params.use_gpu.toString() == 'true' ? 'process_gpu' : null
    label 'error_ignore'

    container "docker://registry.harmonyhpc.io/baldikacti_umass_edu-pub/chienlab_proteinfold/alphafold3:latest"

    input:
    path("input_a3m/*")
    path af3_db
    path af3_model

    output:
    path ("folds/*")                            , emit: folds
    path ("folds/*/*_summary_confidences.json") , emit: summary_json

    script:
    def args = task.ext.args ?: ''
    """
    mkdir folds
    python3 /app/alphafold/run_alphafold.py --norun_data_pipeline --input_dir input_a3m --output_dir folds --db_dir $af3_db --model_dir $af3_model $args
    """
}