#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=base_space_download
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

ENV_FILE=$1
. $ENV_FILE

BASESPACE_PROJECT_ID=$2

mkdir $FASTQ_DIR

module load singularity
singularity exec --pwd /src \
    --no-home \
    --bind $FASTQ_DIR:/src/data/fastq \
    --bind /etc:/etc \
    --env BASESPACE_ACCESS_TOKEN=$BASESPACE_ACCESS_TOKEN \
    --env BASESPACE_API_SERVER=$BASESPACE_API_SERVER \
    $SINGULARITY_IMAGE \
    bs download project -i $BASESPACE_PROJECT_ID -o /src/data/fastq --extension=fastq.gz