#!/bin/bash

#SBATCH --mem=1G
#SBATCH --cpus-per-task=8

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


# chargement des modules nécessaires
module load fastqc/0.11.9
module load bowtie2/2.3.5
module load samtools/1.9
module load htseq/0.11.3

# numéro de l'échantillons à analyser
sample="SRR2960338"
# répertoire de base
base_dir="/shared/projects/form_2021_29/data/rnaseq_tauri"
# répertoire contenant les fichiers du génome de référence
# (séquence et annotations)
genome_dir="${base_dir}/genome"
# nom du fichier contenant les annotations
annotations="${genome_dir}/GCF_000214015.3_version_140606_genomic_DUO2.gff"
# répertoire contenant les fichiers .fastq.gz
fastq_dir="${base_dir}/reads"


echo "=============================================================="
echo "Contrôle qualité - échantillon ${sample}"
echo "=============================================================="
mkdir -p reads_qc
srun fastqc "${fastq_dir}/${sample}.fastq.gz" --outdir reads_qc

echo "=============================================================="
echo "Alignement des reads sur le génome de référence - échantillon ${sample}"
echo "=============================================================="
srun bowtie2 --threads="${SLURM_CPUS_PER_TASK}" -x "${genome_dir}/O_tauri" -U "${fastq_dir}/${sample}.fastq.gz" -S "map/bowtie-${sample}.sam"

echo "=============================================================="
echo "Conversion en binaire, tri et indexation des reads alignés - échantillon ${sample}"
echo "=============================================================="
srun samtools view --threads="${SLURM_CPUS_PER_TASK}" -b "map/bowtie-${sample}.sam" -o "map/bowtie-${sample}.bam"
srun samtools sort --threads="${SLURM_CPUS_PER_TASK}" "map/bowtie-${sample}.bam" -o "map/bowtie-${sample}.sorted.bam"
srun samtools index "map/bowtie-${sample}.sorted.bam"

echo "=============================================================="
echo "Comptage - échantillon ${sample}"
echo "=============================================================="
mkdir -p count
srun htseq-count --stranded=no --type="gene" --idattr="ID" --order=name --format=bam "map/bowtie-${sample}.sorted.bam" "${annotations}" > "count/count-${sample}.txt"

echo "=============================================================="
echo "Nettoyage des fichiers inutiles - échantillon ${sample}"
echo "=============================================================="
rm -f "map/bowtie-${sample}.sam" "map/bowtie-${sample}.bam"
