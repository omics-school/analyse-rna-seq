#!/bin/bash

#SBATCH --mem=1G

# le script va s'arrêter
# - à la première erreur
# - si une variable n'est pas définie
# - si une erreur est recontrée dans un pipe
set -euo pipefail


module load bowtie2/2.3.5


# répertoire contenant les fichiers du génome de référence
# (séquence et annotations)
genome_dir="/shared/projects/uparis_duo_2020/data/genome"
# nom du fichier contenant le génome de référence
genome_ref="${genome_dir}/GCF_000214015.3_version_140606_genomic.fna"


echo "=============================================================="
echo "Indexation du génome de référence"
echo "=============================================================="
srun bowtie2-build "${genome_ref}" "${genome_dir}/O_tauri"
