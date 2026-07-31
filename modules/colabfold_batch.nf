process COLABFOLD_BATCH {
    tag "$accID"
    label 'gpu'
    label 'error_ignore'

    container "docker://ghcr.io/sokrypton/colabfold:1.5.5-cuda12.2.2"

    input:
    tuple val(accID), path(fasta)
    path ("params/*")
    val  numRec

    output:
    tuple val(accID), path ("*")      , emit: results
    path ("*_toprank.json")           , emit: json
    path ("*.png")                    , emit: multiqc

    script:
    def args = task.ext.args ?: ''
    """
    colabfold_batch \\
        ${fasta} \\
        \$PWD \\
        --num-recycle ${numRec} \\
        --data \$PWD \\
        $args
    
    ln -s *scores_rank_001_*.json ${accID}_toprank.json
    """
}