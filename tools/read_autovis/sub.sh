#! /bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=128GB
#SBATCH --time=96:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

export SINGULARITY_SHELL=/bin/bash

for fname in 1684087370 #1678899080 #1675210948 #1666370606 #1678726283  1678734987  1679615321  1689090392
             
do 

echo $fname   
until singularity exec  ~/containers/katcal.sif ./read_autovis_py3.py ${fname}; do echo ${fname} "### error, retry ###";sleep 10; done

done

