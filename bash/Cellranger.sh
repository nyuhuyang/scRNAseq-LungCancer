#!/bin/bash -l

#SBATCH --job-name=cellranger # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 # CPU cores requested per task (change for threaded jobs)
#SBATCH --mem=128G # Memory needed per node (total)
#SBATCH --error=cellranger_%A_%a.err # File to which STDERR will be written, including job ID
#SBATCH --output=cellranger_%A_%a.out # File to which STDOUT will be written, including job ID

#---------------------Variables to be set-------------------------#
echo "Job ID : $JOB_ID"  ${SLURM_ARRAY_TASK_ID}
PROJECT_NAME="scRNAseq-LungCancer"
path=/athena/elementolab/scratch/yah2014/Projects/${PROJECT_NAME}
fastq_path=${path}/data/fastq
transcriptome="/athena/elementolab/scratch/yah2014/Indexed_genome/refdata-cellranger-mm10-3.0.0"
fastq_dir=$(ls ${fastq_path} | tail -n +${SLURM_ARRAY_TASK_ID}| head -1) # Uses job array for each sample in the folder
localcores=64
localmem=128

cd $path/data
# count
# Repeat this command per sample

cellranger_cmd="
cellranger count \
--id="${fastq_dir}" \
--sample="${fastq_dir}" \
--fastqs=$fastq_path/$fastq_dir  \
--transcriptome="${transcriptome}" \
--localcores="${localcores}" \
--localmem="${localmem}" \
--nosecondary \
--chemistry="SC3Pv3"
"
echo -e "\n CMD: $cellranger_cmd \n"
$cellranger_cmd
