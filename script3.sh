# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail

# numéro des échantillons à analyser
# les numéros sont dans une chaîne de caractères et séparés par un espace
samples="10 41 7"
# nom du fichier contenant le génome de référence
genome=GCF_000214015.3_version_140606_genomic.fna
# nom du fichier contenant les annotations
annotations=GCF_000214015.3_version_140606_genomic_DUO2.gff

for sample in $(echo ${samples})
do
    echo "=============================================================="
    echo "Contrôle qualité - échantillon ${sample}"
    echo "=============================================================="
    fastqc HCA-${sample}_R1.fastq.gz

    echo "=============================================================="
    echo "Indexation du génome de référence"
    echo "=============================================================="
    bowtie2-build ${genome} O_tauri

    echo "=============================================================="
    echo "Alignement des reads sur le génome de référence - échantillon ${sample}"
    echo "=============================================================="
    bowtie2 -x O_tauri -U HCA-${sample}_R1.fastq.gz -S bowtie-${sample}.sam 2> bowtie-${sample}.out

    echo "=============================================================="
    echo "Conversion en binaire, tri et indexation des reads alignés - échantillon ${sample}"
    echo "=============================================================="
    samtools view -b bowtie-${sample}.sam > bowtie-${sample}.bam
    samtools sort bowtie-${sample}.bam -o bowtie-${sample}.sorted.bam
    samtools index bowtie-${sample}.sorted.bam

    echo "=============================================================="
    echo "Comptage - échantillon ${sample}"
    echo "=============================================================="
    htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam bowtie-${sample}.sorted.bam ${annotations} > count-${sample}.txt

    echo "=============================================================="
    echo "Nettoyage des fichiers inutiles - échantillon ${sample}"
    echo "=============================================================="
    rm -f bowtie-${sample}.sam bowtie-${sample}.bam
done
