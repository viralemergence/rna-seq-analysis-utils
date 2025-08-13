#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=star_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --array=0-36:1%10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G

ENV_FILE=$1
. $ENV_FILE

module load singularity
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $FASTP_DIR:/src/fastp \
    --bind $GENOMES_DIR:/src/genomes \
    --bind $STAR_COUNTS_DIR:/src/star_counts \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/gene_count_processing/star_orchestration.py \
    -input_dir /src/fastp \
    -genomes_dir /src/genomes \
    -star_counts_dir /src/star_counts \
    -sample_metadata /src/data/sample_metadata_barebones.csv \
    -array_number $SLURM_ARRAY_TASK_ID