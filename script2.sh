# numéro de l'échantillon à analyser
sample=10
# nom du fichier contenant le génome de référence
genome=GCF_000214015.3_version_140606_genomic.fna
# nom du fichier contenant les annotations
annotations=GCF_000214015.3_version_140606_genomic_DUO2.gff

echo "=============================================================="
echo "Contrôle qualité"
echo "=============================================================="
fastqc HCA-${sample}_R1.fastq.gz

echo "=============================================================="
echo "Indexation du génome de référence"
echo "=============================================================="
bowtie2-build ${genome} O_tauri

echo "=============================================================="
echo "Alignement des reads sur le génome de référence"
echo "=============================================================="
bowtie2 -x O_tauri -U HCA-${sample}_R1.fastq.gz -S bowtie-${sample}.sam > bowtie-${sample}.out 2>&1

echo "=============================================================="
echo "Conversion en binaire, tri et indexation des reads alignés"
echo "=============================================================="
samtools view -b bowtie-${sample}.sam > bowtie-${sample}.bam
samtools sort bowtie-${sample}.bam -o bowtie-${sample}.sorted.bam
samtools index bowtie-${sample}.sorted.bam

echo "=============================================================="
echo "Comptage"
echo "=============================================================="
htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam bowtie-${sample}.sorted.bam ${annotations} > count-${sample}.txt

echo "=============================================================="
echo "Nettoyage des fichiers inutiles"
echo "=============================================================="
rm -f bowtie-${sample}.sam bowtie-${sample}.bam
