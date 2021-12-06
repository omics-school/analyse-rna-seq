# Référence de l'échantillon à analyser.
sample=SRR2960338
# Chemin relatif et nom du fichier contenant le génome de référence.
genome=genome/GCF_000214015.3_version_140606.fna
# Chemin relatif et nom du fichier contenant les annotations.
annotations=genome/GCF_000214015.3_version_140606.gff


echo "=============================================================="
echo "Contrôle qualité"
echo "=============================================================="
mkdir -f reads_qc
fastqc reads/${sample}.fastq.gz --outdir reads_qc

echo "=============================================================="
echo "Indexation du génome de référence"
echo "=============================================================="
mkdir -p index
bowtie2-build ${genome} index/O_tauri

echo "=============================================================="
echo "Alignement des reads sur le génome de référence"
echo "=============================================================="
mkdir -p map
bowtie2 -p 2 -x index/O_tauri -U reads/${sample}.fastq.gz -S map/bowtie-${sample}.sam 2> map/bowtie-${sample}.out

echo "=============================================================="
echo "Conversion en binaire, tri et indexation des reads alignés"
echo "=============================================================="
samtools view -@ 2 -b map/bowtie-${sample}.sam -o map/bowtie-${sample}.bam
samtools sort -@ 2 map/bowtie-${sample}.bam -o map/bowtie-${sample}.sorted.bam
samtools index map/bowtie-${sample}.sorted.bam

echo "=============================================================="
echo "Comptage"
echo "=============================================================="
mkdir -p count
htseq-count --stranded=no --type="gene" --idattr="ID" --order=name --format=bam map/bowtie-${sample}.sorted.bam ${annotations} > count/count-${sample}.txt

echo "=============================================================="
echo "Nettoyage des fichiers inutiles"
echo "=============================================================="
rm -f map/bowtie-${sample}.sam map/bowtie-${sample}.bam
