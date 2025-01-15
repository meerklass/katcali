#! /bin/bash

cat result_level1_list.txt | while read line

do

fname=`echo $line | awk '{print $1}'`
ant=`echo $line | awk '{print $2}'`

for pol in v h

do 

echo ${fname} ${ant}  ${pol}

echo "#! /bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=64GB
#SBATCH --time=72:00:00
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out

export SINGULARITY_SHELL=/bin/bash" > KATcali_${fname}_${ant}_${pol}

echo "singularity exec ~/containers/katcal.sif ./KATcali_UHF_level2.py ${fname} ${ant} $pol" >> KATcali_${fname}_${ant}_${pol}

sbatch KATcali_${fname}_${ant}_${pol}

rm -f KATcali_${fname}_${ant}_${pol}
   
done

done
