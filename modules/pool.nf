process POOL {
    label 'process_low'

    container "docker://baldikacti/chienlab_proteinfold_py:latest"

    input:
    tuple path(bait_file), path(pool_file), path(init_pools)
    val export_mode

    output:
    path ("pools/pool_*.{fasta,yaml,json}"), emit: pools
    path ("pools/pools.tsv")               , emit: pools_tsv

    script:
    def args = task.ext.args ?: ''
    def bait_arg = bait_file ? "--bait-fasta ${bait_file}" : ''
    def init_arg = init_pools ? "--init-pools ${init_pools}" : ''
    """
    afpool.py \\
        --pool-fasta ${pool_file} \\
        --output pools \\
        --export-mode ${export_mode} \\
        ${bait_arg} \\
        ${init_arg} \\
        ${args}
    """
}

// Build the POOL input channel from the pool samplesheet.
def poolInput(accession_file) {

    def ch_entries = accession_file
        .splitCsv(header: true, sep: '\t', quote: '"')
        .map { row -> tuple(file(row.Entry, checkIfExists: true), row.bait as Integer) }
        .branch { entry, bait ->
            bait: bait == 1
                return entry
            not_bait: bait == 0
                return entry
        }

    // `--init-pools` is only read in all_vs_all mode
    def init_pools = params.init_pool && params.pool_mode == 'all_vs_all'
        ? file(params.init_pool, checkIfExists: true)
        : []

    return params.pool_mode == 'all_vs_all'
        ? ch_entries.not_bait.map { pool -> tuple([], pool, init_pools) }
        : ch_entries.bait.combine(ch_entries.not_bait).map { bait, pool -> tuple(bait, pool, init_pools) }
}
