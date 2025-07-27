process PREPROCESS_AF3 {
    label 'process_single'
    publishDir "${params.outdir}/json", mode: 'copy', pattern: '*.json'

    container "docker://baldikacti/chienlab_proteinfold_py:latest"

    input:
    path acc_file

    output:
    path ("*.json") , emit: json

    script:
    """
    tsv2json.py --output-dir . --workdir ${workflow.launchDir} ${acc_file}
    """
}