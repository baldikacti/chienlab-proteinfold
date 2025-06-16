# Version v0.6.0

- Allows `fasta` file inputs with one or more entries in the input `tsv` file.

- Allows for custom `--host_url` for ColabFold MSAserver.

- `--org_ref` paramater is not optional.

- Fixed certain edge cases.

# Version v0.5.1

- Fixed a bug in fasta pair generation script that failed to catch download errors.

# Version v0.5.0

- Custom fasta files can be used as bait.
- Fixed a bug that prevented relative paths to fail.
- Added new `-profile`'s for local execution. (Docker, Apptainer, Singularity)

# Version v0.4.0

- GPU resource request optimization for slurm. Now uses better defaults for more optimized performance.

# Version v0.3.0

- Adds a new paramater, `top_rank`, that reruns the top selected number of predictions with 20 recycles and exports them to a new folder called `toprank`.
- Added CHANGELOG.md
- Cleans up the `nextflow.config` to use the `unity` profile on `nf-core/configs`.

# Version v0.2.0

- Fixes a bug that made each colabfold call to download the colabfold cache over and over again.

# Version v0.1.0

- Initial release
- Contains Colabfold pipeline for the Unity cluster
