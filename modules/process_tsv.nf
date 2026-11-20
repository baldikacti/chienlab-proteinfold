process PROCESS_TSV {
    label 'process_single'

    container "docker://registry.harmonyhpc.io/baldikacti_umass_edu-pub/chienlab_proteinfold/chienlab_proteinfold_py:latest"

    input:
    path acc_file
    val mode

    output:
    path ("*.{fasta,json}") , emit: processed_tsv_output

    script:
    """
    tsv2json.py --output-dir . --workdir ${workflow.launchDir} --mode ${mode} ${acc_file}
    """
}