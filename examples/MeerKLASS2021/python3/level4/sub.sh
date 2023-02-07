#! /bin/bash

cat result_level3_list.txt | while read line

do

fname=`echo $line | awk '{print $1}'`
ant=`echo $line | awk '{print $2}'`

echo ${fname} ${ant} 

echo "#! /bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=1:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

export SINGULARITY_SHELL=/bin/bash" > sub_${fname}_${ant}

echo "singularity exec ~/containers/katcal.sif ./KATcali_multi_level4_py3.py ${fname} ${ant}" >> sub_${fname}_${ant}

sbatch sub_${fname}_${ant}

rm -f sub_${fname}_${ant}
   
done
