#!/bin/bash -l
#SBATCH --job-name=hisat2-array
#SBATCH --array=1-192
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=1-05:005:00
#SBATCH --partition brc
#SBATCH --mem=6000MB
#SBATCH --mail-type=END,FAIL
# Autogenerated script from write_hisat2_script.R
# date Fri Jun  5 14:58:54 2020
# make sure directory paths exist before running script



module load apps/hisat2/2.1.0-python3.7.3



# Parse parameter file to get variables.
number=$SLURM_ARRAY_TASK_ID
paramfile=hisat2_param.txt
 
inr1=`sed -n ${number}p $paramfile | awk '{print $1}'`
inr2=`sed -n ${number}p $paramfile | awk '{print $2}'`
outsam=`sed -n ${number}p $paramfile | awk '{print $3}'`

# 9. Run the program.
hisat2 --trim5 10 --threads 4 -x /users/k1625253/scratch/old-scratch_rosalind-legacy-import_2020-01-28/Data/MetaData/GenomeIndex/mm10_hisat2/genome -1 $inr1 -2 $inr2 -S $outsam



