process FASTA2JSON {
    label 'process_single'

    container "docker://baldikacti/chienlab_proteinfold_py:latest"

    input:
    path ("fasta/*")

    output:
    path ("results/*.json"), emit: af3_json

    script:
    """
    fasta2json.py -i fasta -o results
    """
}
