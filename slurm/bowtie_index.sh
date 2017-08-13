#!/bin/bash
#SBATCH --nodes=1               # 1 node
#SBATCH --cpus-per-task=2       # 2CPU cores
#SBATCH --time=1:00:00          # kill after 1 hour
#SBATCH --mem=10GB              # use 10GB of RAM
#SBATCH --job-name=Bowtie_Index # name of the job
#SBATCH --mail-type=ALL         # ALL = Begin, End, and Fail
##SBATCH --mail-user=%u@nyu.edu # this is automatically the submitting
                                # user and is not needed
#SBATCH --output=slurm_%j.out   # put STERR and STDOUT into file called 
                                # slurm_jobid.out

module load bowtie/gnu/1.2.0

python /scratch/nah9/src/dhj_c/dhj_c/bowtie_index.py \
    -f /scratch/nah9/data/S288c_HindIII_adj.fasta \
    -o /scratch/nah9/bowtie_indexes/S288c_HindIII
python /scratch/nah9/src/dhj_c/dhj_c/bowtie_index.py \
    -f /scratch/nah9/data/SK1_HindIII_adj.fasta \
    -o /scratch/nah9/bowtie_indexes/SK1_HindIII

exit 0;



