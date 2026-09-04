# Site-specific deployment

`submit_cyripanel.slurm` is the SLURM submission script used by the clinical
pharmacogenomics service at National Taiwan University Hospital, running on the
NCHC cluster. It hard-codes that site's account, partition and paths.

It is kept here as a worked example of deploying CyriPanel behind a scheduler.
Adapt the `#SBATCH` directives and the conda activation line to your own site;
the tool itself is invoked exactly as in `run_cyripanel.sh` at the repository
root, which is portable.
