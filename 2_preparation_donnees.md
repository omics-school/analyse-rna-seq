---
title: Préparation des données RNA-seq
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

L'objectif de cette partie est de télécharger et contrôler les données RNA-seq de *O. tauri* nécessaires à l'analyse.

## 2.1 Sélectionner les données RNA-seq

Le jeu de données initial, publié en [2016](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2666-6), est constitué de 47 fichiers *.fastq*. Ce projet porte l'identifiant [PRJNA304086](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA304086) sur NCBI / SRA.

Pour que cette activité se déroule dans un temps raisonnable, nous travaillerons sur un jeu de données réduit correspondant à la condition « S2 » (voir la [figure 1](https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-016-2666-6/MediaObjects/12864_2016_2666_MOESM1_ESM.pptx) dans les *additional files*), c'est-à-dire :
- « *Condition 1 Short-term adaptative response of cells* »,
- *Iron (+ Fe)*,
- *Light*,
- *3 hours*.

La [liste des échantillons](sample_info_v2.txt) nous indique que les échantillons associés à ce jeu de données réduit sont *HCA.3*, *HCA.4* et *HCA.5*.

Le site [SRA Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/) est utile pour retrouver les noms des *runs* correspondants.

- Entrez le numéro du projet : PRJNA304086, puis cliquez sur le bouton *Search*
- Repérez les noms des *Run* qui correspondent au *Sample Name* (ou *Library Name*), ici *HCA.3*, *HCA.4* et *HCA.5*.
- Vous devriez trouver les identifiants *SRR2960338*, *SRR2960341* et *SRR2960343*. Notes ces identifiants car vous en aurez besoin par la suite.


## 2.2 Préparer l'arborescence

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

*Explications : cette commande crée les répertoires `rnaseq_tauri`, `rnaseq_tauri/genome` et `rnaseq_tauri/reads` en une seule opération. L'option `-p` ne génère pas d'erreur si les répertoires existent déjà et autorise la création de répertoires imbriqués (par exemple le répertoire `genome` dans le répertoire `rnaseq_tauri`).*

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

Si ce n'est pas le cas, vérifiez que vous n'avez pas fait d'erreur de typo.

## 2.3 Télécharger les fichiers Fastq

Si cela n'est pas déjà fait, activez ensuite l'environnement conda *rnaseq-env* qui contient tous les outils dont vous avez besoin :

```bash
$ conda activate rnaseq-env
```

Téléchargez les 3 fichiers fastq avec la commande suivante :

```bash
$ fasterq-dump --progress --outdir reads SRR2960338 SRR2960341 SRR2960343
```

La même commande aurait pu s'écrire :

```bash
$ fasterq-dump -p -O reads SRR2960338 SRR2960341 SRR2960343
```

L'utilisation des versions longues des options est souvent plus explicite. Essayez de comprendre la signification de ces différentes options. N'hésitez pas à consulter l'aide de `fasterq-dump` en lançant la commande `fasterq-dump --help` dans un autre terminal (après avoir chargé l'environnement conda *rnaseq-env*).

Le téléchargement des données va prendre 20 à 25 minutes. Soyez patient et profitez-en pour prendre un café ou un thé.

Calculez la taille occupée par les fichiers téléchargés avec la commande `du` :

```bash
$ du -ch reads/*
```

*Explications : la commande `du` affiche la taille occupée par des fichiers. L'option `-h` affiche la taille en ko, Mo, Go... L'option `-c` calcule la taille totale occupée par tous les fichiers.*

Les fichiers fastq occupent plus de 5 Go de données ce qui est conséquent pour seulement 3 fichiers. Nous allons les compresser pour gagner un peu de place :

```bash
$ gzip reads/*
```

Cette commande va prendre 5 à 6 minutes. Le temps pour un nouveau café ou thé ?

*Remarque : la commande `gzip` n'affiche rien. C'est normal.*

Vérifiez maintenant le gain obtenu :

```bash
$ du -ch reads/*
```

L'espace disque occupé est désormais inférieur à 1,5 Go pour les 3 fichiers, ce qui est déjà plus raisonnable.

*Remarque : cette étape de compression des données n'est pertinente que parce que les outils que nous utiliserons ensuite savent manipuler ce genre de fichiers.*

## 2.4 Télécharger les fichiers Fastq (solution alternative)

Nous vous présentons ici une solution alternative pour télécharger des fichiers fastq depuis un miroir de SRA situé à l'EBI.

Nous allons utiliser le site [SRA EXplorer](https://sra-explorer.info/) qui est très pratique.

- Sur ce site, indiquez d'abord le numéro du projet, ici PRJNA304086, puis cliquez sur le petite loupe pour lancer la recherche.
- Vous obtenez ensuite 47 réponses qui correspondent au 47 fichiers / échantillons.
- Vous pouvez affiner les réponses en tapant par exemple « *Iron* » dans le champ « *Filter results:* ».
- Sélectionnez les 3 échantillons correspondant à « *Condition 1, Iron, Light, 3H* ».
- Cliquez sur le bouton « *Add 3 to collection* ».
- Cliquez ensuite en haut à droite sur le bouton « *3 saved datasets* ».
- Cliquez enfin sur « *Bash script for downloading FastQ files* ».

Vous obtenez un script Bash qui contient les commandes pour télécharger directement vos fichiers fastq compressés :

```bash
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/003/SRR2960343/SRR2960343.fastq.gz -o SRR2960343_Condition1_Iron_Light_3H.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/001/SRR2960341/SRR2960341.fastq.gz -o SRR2960341_Condition_1_Iron_Light_3H.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/008/SRR2960338/SRR2960338.fastq.gz -o SRR2960338_Condition_1_Iron_Light_3H.fastq.gz
```

⚠️ **Ne lancez pas les commandes ci-dessus car vous avez déjà téléchargé vos données.** ⚠️

Remarque : nous reviendrons sur la notion de script un peu plus tard.

## 2.5 Télécharger le génome de référence

Le génome de référence de *Ostreococcus tauri* et ses annotations sont disponibles sur la [page dédiée](https://www.ncbi.nlm.nih.gov/genome/373?genome_assembly_id=352933) sur le site du NCBI :
- [génome de référence](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.fna.gz)
- [annotations](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.gff.gz). Nous avons légèrement modifié le fichier d'annotations pour ne prendre en compte que les gènes et alléger la visualisation dans IGV.

Vérifiez que vous êtes toujours dans le bon répertoire :

```bash
$ pwd
/mnt/c/Users/omics/rnaseq_tauri
```

Si ce n'est pas le cas, déplacez-vous au bon endroit.

Téléchargez ensuite le génome de référence :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/GCF_000214015.3_version_140606.fna -P genome/
```

L'option `-P` précise où enregistrer le fichier téléchargé.

Les annotations :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/GCF_000214015.3_version_140606.gff -P genome/
```

Et enfin un fichier de sommes de contrôle :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/md5sum.txt -P genome/
```


## 2.6 Contrôler l'intégrité des fichiers

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

L'outil `md5sum` va lire ce fichier, calculer la somme de contrôle du fichier indiqué dans la seconde colonne et la comparer à celle précisée dans la première colonne. Si la somme de contrôle calculée et celle lue correspondent alors `md5sum` affiche « OK » ou « Réussi ».

Normalement, vous devriez obtenir cela :


```bash
$ md5sum -c genome/md5sum.txt
genome/GCF_000214015.3_version_140606.gff: OK
genome/GCF_000214015.3_version_140606.fna: OK
```

Par comparaison avec le contenu du fichier `genome/md5sum.txt`, on peut conclure que l'intégrité des fichiers `genome/GCF_000214015.3_version_140606.fna` et `genome/GCF_000214015.3_version_140606.gff` est vérifiée. Nous avons donc téléchargé les bons fichiers. 🎉

Nous avons vérifié ici l'intégrité de 2 fichiers en une opération, mais il est possible de faire la même chose pour des dizaines de fichiers.

Dans cette étape, vous avez préparé votre arboscence de travail et téléchargé les données nécessaires pour votre analyses RNA-seq. Félicitation 🥳

Vous pouvez fermer votre terminal et passer à l'étape suivante : [➡️](3_analyse_RNA-seq.md)