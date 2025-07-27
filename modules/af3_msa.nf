process AF3_MSA {
    tag "$id"
    label "process_high"
    publishDir "${params.outdir}/${params.mode}/msa/$id", mode: 'copy', pattern: "${id}_data.json"

    input:
    tuple val(id), path(json)
    path af3_db

    output:
    tuple val(id), path("${id}_data.json"), emit: af3_json_processed

    script:
    """
    module load alphafold3/latest
    run_alphafold.py --norun_inference --json_path $json --output_dir msa_results --db_dir $af3_db
    ln -s \$(find msa_results -name "*.json" -type f) ${id}_data.json
    """
}