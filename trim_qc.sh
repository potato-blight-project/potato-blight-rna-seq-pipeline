#!/bin/bash
#SBATCH --job-name trim_qc
#SBATCH --time 240:00:00
#SBATCH --mail-type ALL
#SBATCH --nodes 1
#SBATCH --ntasks 20
#SBATCH --array 1-36
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --account=leachlj-potato-qtl-project

# Load necessary modules
# ----------------------

set -e
module purge; module load bluebear
module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4

# Define the reads directory and sample prefix depending on the sample ID
# -----------------------------------------------------------------------

if [[ "${SLURM_ARRAY_TASK_ID}" -le 18 ]]; then
	READS_DIRECTORY="/rds/projects/l/leachlj-potato-qtl-project/raw_data_blight/S1_18/X201SC21043068-Z01-F001/raw_data"
	SAMPLE_PREFIX="R"
else
	READS_DIRECTORY="/rds/projects/l/leachlj-potato-qtl-project/raw_data_blight/S19_36/X201SC21043068-Z01-F003/raw_data"
	SAMPLE_PREFIX="A"
fi

# Define the output directory and the sample to process
# -----------------------------------------------------

OUTPUT_DIRECTORY="/rds/projects/l/leachlj-potato-qtl-project/M6_Area/OE_Scripts/trimmed"
SAMPLE="${SAMPLE_PREFIX}${SLURM_ARRAY_TASK_ID}"

if [[ -d "${OUTPUT_DIRECTORY}" ]]; then
	mkdir -p "${OUTPUT_DIRECTORY}"
fi

# Run trimming with fastqc
# ------------------------

trim_galore \
	-q 26 \
	--phred33 \
	--fastqc \
	--illumina \
	--length 20 \
	-o "${OUTPUT_DIRECTORY}" \
	--paired "${READS_DIRECTORY}/${SAMPLE}/${SAMPLE}_1.fq.gz" "${READS_DIRECTORY}/${SAMPLE}/${SAMPLE}_2.fq.gz"
