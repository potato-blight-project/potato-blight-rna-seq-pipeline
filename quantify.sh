#!/bin/bash
#SBATCH --job-name quantify
#SBATCH --time 96:00:00
#SBATCH --mail-type ALL
#SBATCH --mem 120G
#SBATCH --array 1-36
#SBATCH --output=slurm-%x_%A_%a.out
#SBATCH --account=leachlj-potato-qtl-project

set -e
module purge; module load bluebear
module load HTSeq/0.11.0-foss-2018b-Python-2.7.15

INPUT_DIRECTORY="/rds/projects/l/leachlj-potato-qtl-project/M6_Area/potato-blight-rna-seq-pipeline/aligned"
OUTPUT_DIRECTORY="/rds/projects/l/leachlj-potato-qtl-project/M6_Area/potato-blight-rna-seq-pipeline/quantified"
DATABASE="/rds/projects/l/leachlj-potato-qtl-project/SharedGenomeFilesV6/1.Spud_downloaded"

if [ ! -d "${OUTPUT_DIRECTORY}" ]; then
	mkdir -p "${OUTPUT_DIRECTORY}"
fi

if [[ "${SLURM_ARRAY_TASK_ID}" -le 18 ]]; then
	SAMPLE_PREFIX="R"
else
	SAMPLE_PREFIX="A"
fi

SAMPLE="${SAMPLE_PREFIX}${SLURM_ARRAY_TASK_ID}"

htseq-count \
	-f bam \
	--stranded=no \
	--mode union \
	--nonunique none \
	"${INPUT_DIRECTORY}/${SAMPLE}_Aligned.sortedByCoord.out.bam" \
	"${DATABASE}/potatoV6_hc_gene_models.gtf" > "${OUTPUT_DIRECTORY}/${SAMPLE}_v6_count.txt"
