/*
 * Batch `files` into groups of `size` whose membership and order are fixed by the sorted names of
 * `inputs`, not by the order the files arrive in. Batched processes then get the same inputs on
 * every run, so -resume keeps their cache. Each batch is emitted as soon as all of its files arrive.
 *
 * Each file is matched to its input by baseName, after removing `suffix` from the file's baseName
 * (e.g. '_data' for AF3_MSA's `pool0_data.json`, from input `pool0.json`). To batch the inputs
 * themselves, pass them as `files` too, with an empty suffix.
 */
def batchByInput(inputs, files, size, suffix) {
    // Input baseName -> batch key; groupKey carries the batch size so groupTuple knows when it is full
    def batch_keys = inputs
        .toSortedList { a, b -> a.name <=> b.name }
        .map { sorted ->
            sorted.collate(size).indexed().collectEntries { i, batch ->
                batch.collectEntries { f -> [f.baseName, groupKey(i, batch.size())] }
            }
        }

    return files
        .combine(batch_keys)
        .map { f, keys -> [keys[f.baseName.replaceFirst(suffix + '$', '')], f] }
        // remainder: a batch missing a file (e.g. from an ignored failure) still runs with the rest
        .groupTuple(remainder: true, sort: { a, b -> a.name <=> b.name })
        .map { _batch, batch_files -> batch_files }
}
