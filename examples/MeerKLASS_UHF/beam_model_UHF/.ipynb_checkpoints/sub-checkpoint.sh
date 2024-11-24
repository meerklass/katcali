#! /bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64000
#SBATCH --time=20:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

export SINGULARITY_SHELL=/bin/bash

singularity exec ~/containers/katcal.sif ./fit_beam_UHF.py


