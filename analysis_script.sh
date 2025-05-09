#!/bin/sh
#SBATCH --job-name=preprocess_reads		# Job name
#SBATCH --mail-type=ALL   				# Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=your_email_here		# Where to send mail
#SBATCH --cpus-per-task=				# Number of cores
#SBATCH --mem-per-cpu=5Gb				# Per processor memory
#SBATCH -t 06:00:00						# Walltime
#SBATCH -o preprocessing.out			# Name output file
#SBATCH --account=barbazuk-b			# Project to run under
#SBATCH --qos=barbazuk-b				# Queue type
#

start_date_time="`date "+%Y-%m-%d %H:%M:%S"`"
echo "Starting pipeline. Started on $start_date_time"

## 1. Quality check & trim reads to remove adapter sequences and low quality reads
echo "Quality checking and trimming reads to remove adapter sequences and low quality reads."

cd /your/directory/here #Make sure all fastq.gz files are placed in this directory

if [! -f ./TruSeq3-PE.fa]; then
	echo "TruSeq3-PE.fa not found. Download TruSeq3-PE.fa before continuing. Cancelling job."
	scancel $SLURM_JOB_ID
fi

module load trimmomatic
module load fastqc

num_files=$(ls *_1P.fastq | wc -l) # Get number of samples
f1=$(ls *_1P.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p) # Get name of each input file
f2=$(ls *_2P.fastq | sed -n ${SLURM_ARRAY_TASK_ID}p)

trimmomatic PE 