process COLABFOLD_BATCH {
    tag "$accID"
    label 'process_medium'
    // toString() guards against a command-line `--use_gpu false`, which reaches a
    // process directive as the String 'false' rather than a boolean
    label params.use_gpu.toString() == 'true' ? 'process_gpu' : null
    label 'error_ignore'

    container "docker://ghcr.io/sokrypton/colabfold:main-cuda13"

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