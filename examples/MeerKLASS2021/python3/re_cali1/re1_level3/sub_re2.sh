#! /bin/bash

cat result_level2_list.txt | while read line

do

fname=`echo $line | awk '{print $1}'`
ant=`echo $line | awk '{print $2}'`
pol=`echo $line | awk '{print $3}'`

echo ${fname} ${ant}  ${pol}

echo "#! /bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=16GB
#SBATCH --time=72:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH -x compute-102

export SINGULARITY_SHELL=/bin/bash" > sub_${fname}_${ant}_${pol}

echo "singularity exec ~/containers/katcal.sif ./KATcali_multi_level3_re2.py ${fname} ${ant} ${pol}" >> sub_${fname}_${ant}_${pol}

sbatch sub_${fname}_${ant}_${pol}

rm -f sub_${fname}_${ant}_${pol}
   
done
