process AF3_MSA {
    label "process_high"
    publishDir "${params.outdir}/${params.mode}", mode: 'copy', pattern: "msa/*.json"

    input:
    path(json)
    path af3_db

    output:
    path("msa/*.json"), emit: af3_json_processed

    script:
    """
    module load alphafold3/latest
    run_alphafold.py --norun_inference --json_path $json --output_dir msa --db_dir $af3_db
    """
}