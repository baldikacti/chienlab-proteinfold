process PROCESS_TSV {
    label 'process_single'
    publishDir "${params.outdir}/${params.mode}/preprocessing", mode: 'copy', pattern: '*.{fasta,json}'

    container "docker://baldikacti/chienlab_proteinfold_py:latest"

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