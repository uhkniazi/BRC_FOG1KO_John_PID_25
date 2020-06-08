#!/bin/bash -l
#SBATCH --job-name=samtools-array
#SBATCH --array=1-192
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-05:05:00
#SBATCH --partition brc
#SBATCH --mem-per-cpu=2000MB
#SBATCH --mail-type=END,FAIL
# Autogenerated script from write_samtools_script.R
# date Mon Jun  8 11:49:02 2020
# make sure directory paths exist before running script



module load apps/samtools/1.9.0-singularity



# Parse parameter file to get variables.
           number=$SLURM_ARRAY_TASK_ID
           paramfile=samtools_param.txt
           
           ins1=`sed -n ${number}p $paramfile | awk '{print $1}'`
           bamfile=`sed -n ${number}p $paramfile | awk '{print $2}'`
           bamq10=`sed -n ${number}p $paramfile | awk '{print $3}'`
           bamq10sort=`sed -n ${number}p $paramfile | awk '{print $4}'`
           bamrd=`sed -n ${number}p $paramfile | awk '{print $5}'`
           
           # 9. Run the program.
samtools view -b -S $ins1 > $bamfile
samtools view -b -q 10 $bamfile > $bamq10
samtools sort -o $bamq10sort $bamq10
samtools index $bamq10sort
samtools rmdup $bamq10sort $bamrd
samtools index $bamrd



