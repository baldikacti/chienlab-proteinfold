process AF3_MSA {
    label 'process_high'
    label 'error_ignore'

    container "docker://registry.harmonyhpc.io/baldikacti_umass_edu-pub/chienlab_proteinfold/alphafold3:latest"

    input:
    path(json)
    path af3_db

    output:
    path("*_data.json"), emit: af3_json_processed

    script:
    def name = json.baseName
    """
    python3 /app/alphafold/run_alphafold.py --norun_inference --json_path $json --output_dir msa --db_dir $af3_db
    ln -s \$(find msa -name "*.json" -type f) ${name}_data.json
    """
}