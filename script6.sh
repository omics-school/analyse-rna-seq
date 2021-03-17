#!/bin/bash

#SBATCH --mem=1G
#SBATCH --cpus-per-task=8
#SBATCH --array=0-46            # 47 samples in total

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
sample="10"
# répertoire contenant les fichiers du génome de référence
# (séquence et annotations)
genome_dir="/shared/projects/uparis_duo_2020/data/genome"
# nom du fichier contenant les annotations
annotations="${genome_dir}/GCF_000214015.3_version_140606_genomic_DUO2.gff"
# répertoire contenant les fichiers .fastq.gz
fastq_dir="/shared/projects/uparis_duo_2020/data/reads"
# listes de tous les fichiers .fastq.gz
fastq_files=(${fastq_dir}/*fastq.gz)
# extraction du numéro de l'échantillon
sample=$(basename -s _R1.fastq.gz "${fastq_files[$SLURM_ARRAY_TASK_ID]}" | cut -d "-" -f 2)

echo "=============================================================="
echo "Contrôle qualité - échantillon ${sample}"
echo "=============================================================="
srun fastqc "${fastq_dir}/HCA-${sample}_R1.fastq.gz" -o .

echo "=============================================================="
echo "Alignement des reads sur le génome de référence - échantillon ${sample}"
echo "=============================================================="
srun bowtie2 --threads="${SLURM_CPUS_PER_TASK}" -x "${genome_dir}/O_tauri" -U "${fastq_dir}/HCA-${sample}_R1.fastq.gz" -S "bowtie-${sample}.sam"

echo "=============================================================="
echo "Conversion en binaire, tri et indexation des reads alignés - échantillon ${sample}"
echo "=============================================================="
srun samtools view -b "bowtie-${sample}.sam" -o "bowtie-${sample}.bam"
srun samtools sort "bowtie-${sample}.bam" -o "bowtie-${sample}.sorted.bam"
srun samtools index "bowtie-${sample}.sorted.bam"

echo "=============================================================="
echo "Comptage - échantillon ${sample}"
echo "=============================================================="
srun htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam "bowtie-${sample}.sorted.bam" "${annotations}" > count-${sample}.txt

echo "=============================================================="
echo "Nettoyage des fichiers inutiles - échantillon ${sample}"
echo "=============================================================="
rm -f "bowtie-${sample}.sam" "bowtie-${sample}.bam" "bowtie-${sample}.sorted.bam.bai"
