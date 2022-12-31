#!/bin/bash

#SBATCH --job-name=GenoPrimer
#SBATCH --time=14-00:00:00
#SBATCH --array=1-39866%1000
#SBATCH --partition preempted
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1
#SBATCH -e slurm.out/slurm-%A_%a.err
#SBATCH -o slurm.out/slurm-%A_%a.out

module load anaconda
conda activate GenoPrimer

#working dir
working_dir="/hpc/projects/data_lg/duo.peng/GenoPrimer/batch"
cd $working_dir

# declare arrays and generate input
readarray -t input < <(cat canonical_N_C_sites_list_GRCh38.csv)
declare -x idx=$(( ${SLURM_ARRAY_TASK_ID}))  #index start from 1  (0=header)

line=${input[$idx]}
IFS=',' read -ra fields <<< "$line"
Genome=${fields[1]}
Chrom=${fields[2]}
Pos=${fields[4]}
Type=${fields[5]}

#main command
cmd="python /hpc/projects/data_lg/duo.peng/GenoPrimer/GenoPrimer.py "
cmd+="--oneliner_input ${Genome},${Chrom},${Pos} "
cmd+="--outdir batch/${Genome}/${Chrom}/${Pos}_${Type} " #output dir is relative to GenoPrimer.py

$cmd
