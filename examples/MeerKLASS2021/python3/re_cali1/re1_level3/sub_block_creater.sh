#! /bin/bash

cat result_level2_list_h.txt | while read line

do

fname=`echo $line | awk '{print $1}'`
ant=`echo $line | awk '{print $2}'`

echo ${fname} ${ant}

echo "#! /bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=2:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH -x compute-102

export SINGULARITY_SHELL=/bin/bash" > sub_${fname}_${ant}

echo "singularity exec ~/containers/katcal.sif ./Tsky_block_creater.py ${fname} ${ant} " >> sub_${fname}_${ant}

sbatch sub_${fname}_${ant}

rm -f sub_${fname}_${ant}
   
done
