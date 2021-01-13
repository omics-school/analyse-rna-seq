# Le script va s'arrêter :
# - à la première erreur,
# - si une variable n'est pas définie,
# - si une erreur est recontrée dans un pipe.
set -euo pipefail

# Numéro des échantillons à analyser.
# Les numéros sont entre guillemets et séparés par un espace.
# Faites en sorte que ces numéros correspondent à VOS échantillons.
samples="10 41 7"
# Chemin relatif et nom du fichier contenant le génome de référence.
genome=genome/GCF_000214015.3_version_140606.fna
# Chemin relatif et nom du fichier contenant les annotations.
annotations=genome/GCF_000214015.3_version_140606.gff


# On indexe le génome qu'une seule fois.
echo "=============================================================="
echo "Indexation du génome de référence"
echo "=============================================================="
mkdir -p index
bowtie2-build ${genome} index/O_tauri


for sample in ${samples}
do
    echo "=============================================================="
    echo "Contrôle qualité - échantillon ${sample}"
    echo "=============================================================="
    fastqc reads/HCA-${sample}_R1.fastq.gz

    echo "=============================================================="
    echo "Alignement des reads sur le génome de référence - échantillon ${sample}"
    echo "=============================================================="
    mkdir -p map
    bowtie2 -p 2 -x index/O_tauri -U reads/HCA-${sample}_R1.fastq.gz -S map/bowtie-${sample}.sam 2> map/bowtie-${sample}.out

    echo "=============================================================="
    echo "Conversion en binaire, tri et indexation des reads alignés - échantillon ${sample}"
    echo "=============================================================="
    samtools view -@ 2 -b map/bowtie-${sample}.sam > map/bowtie-${sample}.bam
    samtools sort -@ 2 map/bowtie-${sample}.bam -o map/bowtie-${sample}.sorted.bam
    samtools index map/bowtie-${sample}.sorted.bam

    echo "=============================================================="
    echo "Comptage - échantillon ${sample}"
    echo "=============================================================="
    mkdir -p count
    htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam map/bowtie-${sample}.sorted.bam ${annotations} > count/count-${sample}.txt

    echo "=============================================================="
    echo "Nettoyage des fichiers inutiles - échantillon ${sample}"
    echo "=============================================================="
    rm -f map/bowtie-${sample}.sam map/bowtie-${sample}.bam
done
