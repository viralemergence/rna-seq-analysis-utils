singularity pull rsau.sif docker://ghcr.io/viralemergence/rna-seq-analysis-utils:latest

export myscratch="$(mkworkspace)"
echo $myscratch

sbatch base_space_download.sh .env 423121844
sbatch base_space_download.sh .env 432190937

# Collate individual basespace fastq.gz into single folder
mkdir $FASTQ_COLLATED_DIR
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $FASTQ_DIR:/src/fastq \
    --bind $FASTQ_COLLATED_DIR:/src/fastq_collated \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/read_preprocessing/basespace_read_collation.py \
    -basespace_dir /src/fastq \
    -collated_dir /src/fastq_collated

# Run fastp on raw sequencing data
mkdir $FASTP_DIR
sbatch fastp_array.sh .env