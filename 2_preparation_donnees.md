---
title: Préparation des données RNA-seq
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

L'objectif de cette partie est de télécharger et contrôler les données RNA-seq de *O. tauri* nécessaires à l'analyse.

## 1 Sélection des données RNA-seq

Le jeu de données initial, publié en [2016](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2666-6), est constitué de 47 fichiers *.fastq.gz* (format *.fastq* compressé) pour un total de 24 Go. Ce projet porte l'identifiant [PRJNA304086](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA304086).

Pour que cette activité se déroule dans un temps raisonnable, nous travaillerons sur un jeu de données réduit correspondant à la condition « S2 », c'est-à-dire :
- « *Condition 1 Short-term adaptative response of cells* »
- *Iron (+ Fe)*
- *Light*
- *3 hours*

Ce jeu de données est composé de 3 fichiers *.fastq.gz* qui correspondent aux 3 réplicats de la condition S2 et qui portent les identifiants :
- SRR2960338
- SRR2960341
- SRR2960343


### 1.1 Préparation de l'arborescence

Sous Windows, ouvrez un terminal Ubuntu. Si vous avez oublié comment faire, consultez le [tutoriel Unix](https://omics-school.github.io/unix-tutorial/tutoriel/README).

Déplacez vous ensuite dans le répertoire `/mnt/c/Users/omics` :

```bash
$ cd /mnt/c/Users/omics
```

🔔 Rappels : 

- Ne tapez pas le `$` en début de ligne et faites attention aux majuscules et aux minuscules (surtout pour `Users`) !
- Utilisez le copier / coller.
- Utilisez la complétion des noms de fichier et de répertoires avec la touche <kbd>Tab</kbd>.

Créez ensuite vos répertoires de travail avec la commande :

```bash
$ mkdir -p rnaseq_tauri/{genome,reads}
```

Déplacez-vous dans le répertoire `rnaseq` :

```bash
$ cd rnaseq_tauri
```

La commande `pwd` doit vous renvoyer :

```bash
$ pwd
/mnt/c/Users/omics/rnaseq_tauri
```

et la commande `tree` :

```bash
$ tree
.
├── genome
└── reads
```

### 1.2 Fichiers Fastq

Activez ensuite l'environnement conda *rnaseq-env* qui contient tous les outils dont vous avez besoin :

```bash
$ conda activate rnaseq-env
```

Téléchargez les 3 fichiers fastq avec la commande suivante :

```bash
$ fasterq-dump --threads 3 --progress --outdir reads SRR2960338 SRR2960341 SRR2960343
```

Le téléchargement des données va prendre plusieurs minutes. Soyez patient et profitez-en pour prendre un café ou un thé.

Calculez la taille occupée par les fichiers de données avec la commande `du` :

```bash
$ du -ch reads/*
```

Explications : la commande `du` affiche la taille occupée par des fichiers. L'option `-h` affiche la taille en ko, Mo, Go... L'option `-c` calcule la taille totale occupée par tous les fichiers.

Les fichiers .fastq occupent plus de 5 Go de données ce qui est assez conséquent. Nous allons les compresser pour gagner un peu de place :

```bash
$ gzip reads/*
```

Cette commande va prendre quelques minutes et faire chauffer votre machine. 
C'est tout à fait normal, les outils de compression sont toujours consommateurs de CPU.

Vérifiez maintenant le gain obtenu :

```bash
$ du -ch reads/*
```

L'espace disque occupé est désormais oins de 1,5 Go pour les 3 fichiers, ce qui est déjà plus raisonnable.

Remarque : cette étape de compression des données n'est pertinente que parce que les outils que nous utiliserons ensuite savent manipuler ce genre de fichiers.

### 1.3 Fichiers Fastq (solution alternative)

Nous vous présentons ici une solution alternative pour télécharger des fichiers fastq depuis SRA.

Ne lancez pas les commandes suivantes car vous avez déjà téléchargé vos données.

Le site [SRA EXplorer](https://sra-explorer.info/) est très pratique.

- Sur ce site, indiquez d'abord le numéro de projet, ici PRJNA304086, puis cliquez sur le petite loupe pour lancer la recherche.
- Vous obtenez ensuite 47 réponses correspond au 47 fichiers.
- Vous pouvez raffiner les réponses en tapant par exemple « *Iron* ».
- Sélectionnez ensuite les 3 échantillons correspondant à « *Condition 1, Iron, Light, 3H* ».
- Cliquez ensuite sur le bouton « *Add 3 to collection* ».
- Cliquez ensuite en haut à droite sur le bouton « *3 saved datasets* »
- Cliquez ensuite sur « *Bash script for downloading FastQ files* ». Vous obtenez un script Bash qui contient les commandes pour télécharger directement vos fichiers fastq compressés.


# 2. Préparation du génome de référence

Le génome de référence de *Ostreococcus tauri* et ses annotations sont disponibles sur la [page dédiée](https://www.ncbi.nlm.nih.gov/genome/373?genome_assembly_id=352933) sur le site du NCBI :
- [génome de référence](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.fna.gz)
- [annotations](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.gff.gz). Nous avons légèrement modifié le fichier d'annotations pour ne prendre en compte que les gènes et alléger la visualisation dans IGV.

## 2.1 Téléchargement des données

Vérifiez que vous ếtes toujours dans le bon répertoire :

```bash
$ pwd
/mnt/c/Users/omics/rnaseq_tauri
```

Si ce n'est pas le cas, déplacez-vous au bon endroit.

Téléchargez ensuite le génome de référence :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/GCF_000214015.3_version_140606.fna -P genome/
```

les annotations :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/GCF_000214015.3_version_140606.gff -P genome/
```

et enfin un fichier de somme de contrôle :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/md5sum.txt -P genome/
```

## 2.1 Contrôle de l'intégrité

Affichez le contenu du fichier `genome/md5sum.txt` que vous venez de télécharger avec la commande `cat` :

```bash
$ cat genome/md5sum.txt 
bfef14688f9cbcca45d794ec0348aa2e  genome/GCF_000214015.3_version_140606.gff
046b6e933274c884428a7d5929090f5d  genome/GCF_000214015.3_version_140606.fna
```

Ce fichier contient deux colonnes. Une première avec une somme de contrôle MD5 et une seconde avec un nom de fichier (et son chemin relatif). Il nous permettra de vérifier l'intégrité du fichier contenant le génome de référence et de celui contenant ses annotations en une seule opération :

```bash
$ md5sum -c genome/md5sum.txt
```

L'outil `md5sum` va lire ce fichier, calculer la somme des contrôles des fichiers indiqués dans la seconde colonne et la comparer à celles indquées dans la première colonne. Si la somme de contrôle calculée et celle lue correspondent alors `md5sum` affiche « OK » ou « Réussi ».

Normalement, vous devriez obtenir cela :


```bash
$ md5sum -c genome/md5sum.txt 
genome/GCF_000214015.3_version_140606.gff: Réussi
genome/GCF_000214015.3_version_140606.fna: Réussi
```

Par comparaison avec le contenu du fichier `genome/md5sum.txt`, on peut conclure que l'intégrité des fichiers `genome/GCF_000214015.3_version_140606.fna` et `genome/GCF_000214015.3_version_140606.gff` est vérifiée. Nous avons donc téléchargé les bons fichiers. 🎉

Nous avons vérifié ici l'intégrité de 2 fichiers en une opération, mais il est possible de faire la même chose pour des dizaines de fichiers.

Vous pouvez passer à l'étape suivante : [➡️](3_analyse_RNA-seq.md)