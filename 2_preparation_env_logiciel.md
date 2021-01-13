---
title: Préparation de l'environnement logiciel
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

L'objectif de cette partie est d'installer conda et les logiciels nécessaires à l'analyse RNA-seq des données *O. tauri*.

# Installation de conda

Conda est un gestionnaire de logiciels et d'environnements très utilisé en bioinformatique. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) est une distribution qui permet d'installer conda. Puisque nous travaillons sous Linux (sous Windows 10 certes, mais nous vous rappelons que WSL est un Linux), nous allons installer la version Linux de Miniconda.

Pour cette étape, déplacez-vous dans votre répertoire utilisateur Unix :
```
$ cd 
```

Vérifiez que la commande `pwd` renvoie bien `/home/duo`.

Téléchargez ensuite miniconda :
```
$ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Installez miniconda :
```
$ bash Miniconda3-latest-Linux-x86_64.sh -b -f
```

Puis initialisez conda :
```
$ ./miniconda3/bin/conda init
```

Attention, la commande débute bien par un point `.`

Fermez enfin votre terminal.

Relancez un nouveau terminal. Vous devriez voir `(base)` à gauche de votre invite de commande.

Vérifiez enfin que conda est fonctionnel en tapant la commande :
```
$ conda --version
```

Si vous obtenez `conda 4.9.2` : bravo :tada:


# Création de l'environnement rnaseq

Nous souhaitons maintenant installer les logiciels nécessaires à l'analyse RNA-seq. Nous pourrions le faire dans l'environnement par défaut de Miniconda (qui s'appelle *base* comme l'indique le `(base)` à gauche de votre invite de commande) mais ce serait une très mauvaise pratique.

Nous allons donc créer un environnement conda dédié pour notre analyse RNA-seq. Sans grande originalité, nous allons appelé cet environnement `rnaseq` :

```
$ conda create -n rnaseq -y
```

Une fois l'environnement créé, il faut l'activer (c'est-à-dire l'utiliser) :
```
$ conda activate rnaseq
```

Le `(base)` à gauche de votre invite de commande est maintenant remplacé par `(rnaseq)`.

Voici la liste des logiciels dont nous avons besoin :

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) pour le contrôle qualité des données de séquençage.
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) pour l'indexation du génome de référence puis l'alignement des *reads* sur le génome de référence.
- [SAMtools](http://samtools.sourceforge.net/) pour la manipulation des fichiers d'alignements (conversion en binaire, tri et indexation).
- [HTSeq](https://htseq.readthedocs.io/en/latest/) pour le comptage du nombre de *reads* alignés sur chaque gène.

Nous allons installer ces logiciels dans l'environnement conda *rnaseq* que nous venons de créer :
```
$ conda install -c bioconda fastqc bowtie2 samtools htseq -y
```

L'installation des logiciels et de leurs dépendances va prendre quelques minutes.

Pour terminer, vérifiez que les logiciels sont bien installés en affichant leur version :

```
$ fastqc --version
FastQC v0.11.9
```

```
$ bowtie2 --version
/home/duo/miniconda3/envs/rnaseq/bin/bowtie2-align-s version 2.3.5.1
64-bit
Built on
Wed Apr 17 02:40:25 UTC 2019
[...]
```

```
$ samtools --version
samtools 1.9
Using htslib 1.9
Copyright (C) 2018 Genome Research Ltd.
```

```
$ htseq-count -h
[...]
Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology
Laboratory (EMBL) and Fabio Zanini (fabio.zanini@stanford.edu), Stanford
University. (c) 2010-2019. Released under the terms of the GNU General Public
License v3. Part of the 'HTSeq' framework, version 0.11.3.
```

Bravo :tada: 

Vous avez installé Miniconda, créé un environnement conda et installé tous les logiciels nécessaires.

Vous pouvez passer à l'étape suivante.

