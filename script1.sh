# Numéro de l'échantillon à analyser.
sample=3
# Chemin relatif et nom du fichier contenant le génome de référence.
genome=genome/GCF_000214015.3_version_140606.fna
# Chemin relatif et nom du fichier contenant les annotations.
annotations=genome/GCF_000214015.3_version_140606.gff


echo "Contrôle qualité"
fastqc reads/HCA-${sample}_R1.fastq.gz

echo "Indexation du génome de référence"
mkdir -p index
bowtie2-build ${genome} index/O_tauri

echo "Alignement des reads sur le génome de référence"
mkdir -p map
bowtie2 -p 2 -x index/O_tauri -U reads/HCA-${sample}_R1.fastq.gz -S map/bowtie-${sample}.sam

echo "Conversion en binaire, tri et indexation des reads alignés"
samtools view -@ 2 -b map/bowtie-${sample}.sam > map/bowtie-${sample}.bam
samtools sort -@ 2 map/bowtie-${sample}.bam -o map/bowtie-${sample}.sorted.bam
samtools index map/bowtie-${sample}.sorted.bam

echo "Comptage"
mkdir -p count
htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam map/bowtie-${sample}.sorted.bam ${annotations} > count/count-${sample}.txt

echo "Nettoyage des fichiers inutiles"
rm -f map/bowtie-${sample}.sam map/bowtie-${sample}.bam
