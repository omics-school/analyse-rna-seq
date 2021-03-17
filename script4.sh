#!/bin/bash

#SBATCH --mail-user=your-email@mail.com
#SBATCH --mail-type=ALL
#SBATCH --mem=1G

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


module load fastqc/0.11.9
module load bowtie2/2.3.5
module load samtools/1.9
module load htseq/0.11.3

# numéro des échantillons à analyser
# les numéros sont entre guillemets et séparés par un espace
# faites en sorte que ces numéros correspondant à VOS échantillons
samples="10"
# répertoire contenant les fichiers du génome de référence
# (séquence et annotations)
genome_dir="/shared/projects/uparis_duo_2020/data/genome"
# nom du fichier contenant les annotations
annotations="${genome_dir}/GCF_000214015.3_version_140606_genomic_DUO2.gff"
# répertoire contenant les fichiers .fastq.gz
fastq_dir="/shared/projects/uparis_duo_2020/data/reads"


echo "=============================================================="
echo "Contrôle qualité - échantillon ${sample}"
echo "=============================================================="
srun fastqc "${fastq_dir}/HCA-${sample}_R1.fastq.gz"

echo "=============================================================="
echo "Alignement des reads sur le génome de référence - échantillon ${sample}"
echo "=============================================================="
srun bowtie2 -x "${genome_dir}/O_tauri" -U "${fastq_dir}/HCA-${sample}_R1.fastq.gz" -S "bowtie-${sample}.sam"

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
rm -f bowtie-${sample}.sam bowtie-${sample}.bam
