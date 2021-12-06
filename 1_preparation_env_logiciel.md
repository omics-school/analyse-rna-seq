---
title: Préparation de l'environnement logiciel
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

L'objectif de cette partie est d'installer conda et les logiciels nécessaires à l'analyse RNA-seq des données *O. tauri*.

## 2.1 Installation de conda

Conda est un gestionnaire de logiciels et d'environnements très utilisé en bioinformatique. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) est une distribution qui permet d'installer conda. Puisque nous travaillons sous Linux (sous Windows 10 certes, mais nous vous rappelons que WSL est un système Linux), nous allons installer la version **Linux** de Miniconda.

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

Relancez un nouveau terminal. Vous devriez voir `(base)` à gauche de votre invite de commande (comme [ici](img/conda_base.png)).

Vérifiez enfin que conda est fonctionnel en tapant la commande :

```bash
$ conda --version
```

Si vous obtenez `conda 4.10.3` ou une version supérieure : bravo 🎉


## 2.2 Installation de mamba

Conda est parfois lent à installer un environnement, c'est-à-dire l'ensemble des outils nécessaires pour une tâche particulière, ici une analyse RNA-seq.

Nous vous conseillons d'installer [mamba](https://github.com/mamba-org/mamba) qui vous permettra d'accélérer conda avec la commande :

```bash
$ conda install mamba -n base -c conda-forge
```

## 2.3 Création de l'environnement rnaseq

Nous souhaitons maintenant installer tous les logiciels nécessaires à l'analyse RNA-seq. Nous pourrions le faire dans l'environnement par défaut de Miniconda (qui s'appelle *base* comme l'indique le `(base)` à gauche de votre invite de commande) mais ce serait une très mauvaise pratique.

Nous allons donc créer un environnement conda dédié pour notre analyse RNA-seq. Sans grande originalité, nous appellerons cet environnement `rnaseq` :

```
$ conda create -n rnaseq-env -y
```

Une fois l'environnement créé, il faut l'activer (c'est-à-dire l'utiliser) :
```
$ conda activate rnaseq-env
```

Le `(base)` à gauche de votre invite de commande est maintenant remplacé par `(rnaseq-env)` (comme [ici](img/conda_rnaseq.png)).


## 2.4 Installation des logiciels nécessaires

Voici la liste des logiciels dont nous avons besoin :

- [SRA Toolkit](https://github.com/ncbi/sra-tools) pour télécharger les données de séquençage.
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) pour le contrôle qualité des données de séquençage.
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) pour l'indexation du génome de référence puis l'alignement des *reads* sur le génome de référence.
- [SAMtools](http://samtools.sourceforge.net/) pour la manipulation des fichiers d'alignements (conversion en binaire, tri et indexation).
- [HTSeq](https://htseq.readthedocs.io/en/latest/) pour le comptage du nombre de *reads* alignés sur chaque gène.

Installons ces logiciels dans l'environnement conda *rnaseq* que nous venons de créer. Pour accéler l'installation, nous utilisons ici mamba :
```
$ mamba install -c conda-forge -c bioconda sra-tools fastqc bowtie2 samtools htseq -y
```

L'installation des logiciels et de leurs dépendances va prendre quelques minutes.

Pour terminer, vérifiez que les logiciels sont bien installés en affichant leurs versions :

*Remarque : les versions peuvent légèrement différer.*

```bash
$ fasterq-dump --version

"fasterq-dump" version 2.11.0
```

*Remarque : `fasterq-dump` est un outil fourni par SRA Tools*.


```
$ fastqc --version
FastQC v0.11.9
```

```
$ bowtie2 --version
/home/pierre/.soft/miniconda3/envs/rnaseq/bin/bowtie2-align-s version 2.4.4
64-bit
Built on default-51d38ee7-0be3-4ea2-a3d8-c7558b235b3c
Mon May 24 01:26:39 UTC 2021
[...]
```

```
$ samtools --version
samtools 1.14
Using htslib 1.14
Copyright (C) 2021 Genome Research Ltd.
[...]
```

```
$ htseq-count -h
[...]
Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology Laboratory (EMBL) 
and Fabio Zanini (fabio.zanini@unsw.edu.au), UNSW Sydney. (c) 2010-2020. Released under the terms of
the GNU General Public License v3. Part of the 'HTSeq' framework, version 1.99.2.
```

Bravo ✨

Vous avez installé Miniconda, créé un environnement conda et installé tous les logiciels nécessaires.

Vous pouvez passer à l'étape suivante : [➡️](2_preparation_donnees.md)

