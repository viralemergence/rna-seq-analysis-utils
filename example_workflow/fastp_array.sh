#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=fastp_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --array=0-36:1%10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G

ENV_FILE=$1
. $ENV_FILE

module load singularity
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    --bind $FASTQ_COLLATED_DIR:/src/fastq_collated \
    --bind $FASTP_DIR:/src/fastp \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/read_preprocessing/fastp_orchestration.py \
    -input_dir /src/fastq_collated \
    -output_dir /src/fastp \
    -sample_metadata /src/data/sample_metadata_barebones.csv \
    -array_number $SLURM_ARRAY_TASK_ID \
    -r1_adapter $R1_ADAPTER_SEQ -r2_adapter $R2_ADAPTER_SEQ