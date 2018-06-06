# numéro de l'échantillon à analyser
sample=10
# nom du fichier contenant le génome de référence
genome=GCF_000214015.3_version_140606_genomic.fna
# nom du fichier contenant les annotations
annotations=GCF_000214015.3_version_140606_genomic_DUO2.gff


echo "Contrôle qualité"
fastqc HCA-${sample}_R1.fastq.gz

echo "Indexation du génome de référence"
bowtie2-build ${genome} O_tauri

echo "Alignement des reads sur le génome de référence"
bowtie2 -x O_tauri -U HCA-${sample}_R1.fastq.gz -S bowtie-${sample}.sam

echo "Conversion en binaire, tri et indexation des reads alignés"
samtools view -b bowtie-${sample}.sam > bowtie-${sample}.bam
samtools sort bowtie-${sample}.bam -o bowtie-${sample}.sorted.bam
samtools index bowtie-${sample}.sorted.bam

echo "Comptage"
htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam bowtie-${sample}.sorted.bam ${annotations} > count-${sample}.txt

echo "Nettoyage des fichiers inutiles"
rm -f bowtie-${sample}.sam bowtie-${sample}.bam
