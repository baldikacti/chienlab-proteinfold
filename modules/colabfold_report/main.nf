process COLAB_REPORT {
    label 'process_single'
    publishDir "${params.outdir}/ranked_results", mode: 'copy', pattern: '*.html'

    container "docker://baldikacti/chienlab_proteinfold_rverse:4.4.2"

    input:
    path qmd
    path pngs

    output:
    path ("colabfold_report.html")

    script:
    """
    quarto render ${qmd}
    """
}