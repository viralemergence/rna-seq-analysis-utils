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

# Downloading Rousettus data
mkdir $GENOMES_DIR
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/176/215/GCF_014176215.1_mRouAeg1.p/GCF_014176215.1_mRouAeg1.p_genomic.fna.gz \
    -P $GENOMES_DIR \
    && gzip -d $GENOMES_DIR/GCF_014176215.1_mRouAeg1.p_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/176/215/GCF_014176215.1_mRouAeg1.p/GCF_014176215.1_mRouAeg1.p_genomic.gtf.gz \
    -P $GENOMES_DIR \
    && gzip -d $GENOMES_DIR/GCF_014176215.1_mRouAeg1.p_genomic.gtf.gz

# To generate STAR genome index
singularity exec \
    --pwd /src \
    --no-home \
    --bind $GENOMES_DIR:/src/genomes \
    $SINGULARITY_IMAGE \
    STAR \
    --runThreadN 10 --runMode genomeGenerate \
    --genomeDir /src/genomes/rousettus/ \
    --genomeFastaFiles /src/genomes/GCF_014176215.1_mRouAeg1.p_genomic.fna \
    --sjdbGTFfile /src/genomes/GCF_014176215.1_mRouAeg1.p_genomic.gtf --sjdbOverhang 100 \
    --outFileNamePrefix /src/genomes/rousettus/

# Remove original files
rm $GENOMES_DIR/GCF_014176215.1_mRouAeg1.p_genomic.fna
rm $GENOMES_DIR/GCF_014176215.1_mRouAeg1.p_genomic.gtf

# Get gene counts with STAR
mkdir $STAR_COUNTS_DIR
sbatch star_array.sh .env

# Collate and process gene counts
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $STAR_COUNTS_DIR:/src/star_counts \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/gene_count_processing/gene_count_processing.py \
    -input_dir /src/star_counts \
    -collated_gene_counts /src/data/collated_gene_counts.csv \
    -sample_metadata /src/data/sample_metadata_barebones.csv

# Manually examine counts to determine if anything seems off
# Samples to exclude should be marked in the metadata file by putting 1 in the blacklist column

# Rectify potential batch effect from library preparation
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/gene_count_processing/batch_effect_rectification.py \
    -gene_counts /src/data/collated_gene_counts.csv \
    -sample_metadata /src/data/sample_metadata_barebones.csv \
    -outpath /src/data/gene_counts_final.csv

# Perform correlations between samples per cell line
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/gene_count_correlation_matrix.py \
    -gene_counts /src/data/gene_counts_final.csv \
    -sample_metadata /src/data/sample_metadata_barebones.csv \
    -cell_lines R06E \
    -outdir /src/data/correlations

# Calculate contrasts between viruses, per cell line, per time point
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/calculate_virus_contrasts_per_time.py \
    -gene_counts /src/data/gene_counts_final.csv \
    -sample_metadata /src/data/sample_metadata_barebones.csv \
    -cell_lines R06E \
    -virus_contrasts MR766_vs_No_Virus PRVABC59_vs_No_Virus MR766_vs_PRVABC59 \
    -outdir /src/data/contrasts

# Get "background genes" for gene ontology enrichment analysis (GOEA)
    # Go to: https://www.ncbi.nlm.nih.gov/gene
    # For Rousettus search: "9407"[Taxonomy ID] AND alive[property] AND genetype protein coding[Properties]
    # Send to file (must be tabular)

# Perform coarse GOEA using all time points for a given contrast
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/goatools_analysis_coarse.py \
    -gene_contrast_path /src/data/contrasts/R06E_virus_contrasts_per_time.csv \
    -contrast MR766_vs_No_Virus PRVABC59_vs_No_Virus MR766_vs_PRVABC59 \
    -expression_direction up \
    -taxon_id 9407 \
    -background_genes /src/data/goea/reference_files/ncbi_gene_results_9407.txt \
    -goea_dir /src/data/goea

# Cluster DEGs into patterns
singularity exec --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/deg_patterns/degpattern_orchestration.py \
    -gene_counts /src/data/gene_counts_final.csv \
    -sample_metadata /src/data/sample_metadata_barebones.csv \
    -gene_contrast_path /src/data/contrasts/R06E_virus_contrasts_per_time.csv \
    -cell_line R06E -outdir /src/data/deg_patterns/R_output

# Generate relative gene abundance graphs for gene clusters
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/gene_relative_abundance.py \
    -gene_counts /src/data/gene_counts_final.csv \
    -sample_metadata /src/data/sample_metadata_barebones.csv \
    -gene_clusters_dir /src/data/deg_patterns/R_output \
    -outdir /src/data/deg_patterns

# Perform GOEA on gene clusters
singularity exec \
    --pwd /src \
    --no-home \
    --bind $APP_DIR:/src/app \
    --bind $DATA_DIR:/src/data \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/goatools_analysis_degpatterns.py \
    -gene_clusters_dir /src/data/deg_patterns/R_output \
    -taxon_id 9407 \
    -background_genes /src/data/goea/reference_files/ncbi_gene_results_9407.txt \
    -goea_dir /src/data/goea