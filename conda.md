---
title: Conda et la création d'un environnement de travail pour des analyses reproductibles
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

# Conda et la création d'un environnement de travail pour des analyses reproductibles


## Anaconda, Miniconda, Conda & Bioconda

[Anaconda](https://www.anaconda.com/what-is-anaconda/) est une distribution open source, disponible pour Windows, Mac et Linux, et qui contient de nombreux outils utilisés pour l'analyse de données avec le langage de programmation Python.

Anaconda est basé sur [conda](https://conda.io/docs/), un gestionnaire de paquets, qui permet d'installer des logiciels facilement et sans être administrateur. Conda permet d'installer des logiciels écrits en Python mais aussi en R, en C...

Enfin, Anaconda est également disponible dans une version *light* appelée [Miniconda](https://conda.io/miniconda.html). Miniconda ne contient pas tous les outils Python disponibles dans Anaconda, mais il contient néanmoins le gestionnaire de paquets conda.

Enfin, [Bioconda](https://bioconda.github.io/) est un canal de diffusion de logiciels, utilisable par le gestionnaire de paquets conda et proposant de nombreux logiciels utilisés en bioinformatique.

Voici deux articles très intéressants sur conda :
 
- [Conda le meilleur ami du bioinformaticien](https://bioinfo-fr.net/conda-le-meilleur-ami-du-bioinformaticien). Article d'introduction. Attention cependant, certaines commandes sont obsolètes.
- [Comment fixer les problèmes de déploiement et de durabilité des outils en bioinformatique ? Indice : conda !](https://bioinfo-fr.net/comment-fixer-les-problemes-de-deploiement-et-de-durabilite-des-outils-en-bioinformatique). Article un peu plus technique.

et le papier de référence de Bioconda :

- [Bioconda: sustainable and comprehensive software distribution for the life sciences](https://www.nature.com/articles/s41592-018-0046-7), Björn Grüning et *al.*, Nature methods, 2018.


**Conda est devenu un standard pour installer et utiliser des logiciels en génomique.**


## Installer miniconda dans votre répertoire utilisateur

[Miniconda](https://conda.io/miniconda.html) est une distribution qui contient le gestionnaire de paquets conda.

Téléchargement de la dernière version :

```
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Installation de miniconda :
```
$ bash Miniconda3-latest-Linux-x86_64.sh -b -f
```

Conda sera installé par défaut dans le répertoire `miniconda3` dans votre répertoire utilisateur.

Activation de conda :
```
$ ./miniconda3/bin/conda init
```


Après une déconnexion/reconnexion, on vérifie que conda est bien disponible dans le shell de l'utilisateur avec la commande :
```
$ conda --version
```

## Mettre à jour conda

Il est pertinent de mettre à jour conda après l'installation : 

```
$ conda update -y conda
```

## Gérer un environnement virtuel

Création :
```
$ conda create -n rnaseq
```

Activation de l'environnement :
```
$ conda activate rnaseq
```

Lorsqu'un environnement virtuel est activé, l'invite de commande du shell Linux est modifié et ressemble à :
```
(rnaseq) ppoulain@candihub:~$
```

Pour information, pour quitter un environnement virtuel, il faut utiliser la commande :
```
$ conda deactivate
```


## Installer les logiciels utilisés pour l'analyse RNA-seq :

Il faut avoir pris soin d'activer l'environnement virtuel au préalable.

```
$ conda install -y -c bioconda fastqc bowtie2 htseq samtools
```

*Simple as that!* 😊

L'installation peut parfois prendre un peu de temps.
