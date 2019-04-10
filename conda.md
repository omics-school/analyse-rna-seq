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

## Installation de miniconda

[Miniconda](https://conda.io/miniconda.html) est une distribution qui contient le gestionnaire de paquets conda.

Téléchargement de la dernière version :

```
$ wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Installation de miniconda
```
$ bash Miniconda3-latest-Linux-x86_64.sh -b -p /data/omics-school/share/miniconda
```

`/data/omics-school/share/miniconda` est le répertoire où sera installé miniconda.

Une fois l'installation terminée, il faut ajouter la ligne ci-dessous dans le fichier `.bashrc` qui est situé dans le répertoire utilisateur :
```
source /data/omics-school/share/miniconda/etc/profile.d/conda.sh
```

Rappel : `/data/omics-school/share/miniconda` est le répertoire dans lequel est installé conda.

On vérifie que conda est bien disponible dans le shell de l'utilisateur avec la commande :
```
$ conda --version
```

## Mise-à-jour de conda
```
$ conda update -y conda
```

## Configuration du canal Bioconda

[Bioconda](https://bioconda.github.io/) est un canal de distribution de paquets installables par conda. Bioconda propose de nombreux logiciels dédiés à la bioinfo. La liste complète est disponible [ici](https://bioconda.github.io/recipes.html)

Ajout des canaux `conda-forge` et `bioconda` :
```
$ conda config --add channels conda-forge
$ conda config --add channels bioconda
```

## Gestion de l'environnement virtuel

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

## Installation des logiciels utilisés pour l'analyse RNA-seq :

Il faut avoir pris soin d'activer l'environnement virtuel au préalable.

```
$ conda install -y fastqc bowtie2 htseq samtools
```

*Simple as that!* :smile:


Vérification des versions des logiciels :
```
$ fastqc --version
$ bowtie2 --version | head -n 1
$ samtools --version | head -n 1
$ htseq-count -h | grep version
```


## Remarque pour Ubuntu Server 16.04

Sur Ubuntu Serveur 106.04, par défaut, le fichier `.bashrc` dans le répertoire personnel des utilisateurs n'est pas lu (voir [.bashrc not executed when opening new terminal](https://askubuntu.com/questions/161249/bashrc-not-executed-when-opening-new-terminal))

La solution est alors de créer, pour chaque utilisateur, le fichier `.bash_profile` avec :
```
# include .bashrc if it exists
if [ -f "$HOME/.bashrc" ]; then
    . "$HOME/.bashrc"
fi
```

## Création d'un environnement pour un utilisateur différent de celui qui a installé conda

Conda a été installé dans le répertoire `/data/omics-school/share/miniconda/` par l'utilisateur `ppoulain`.

L'environnement `rnaseq` a été créé par le même utilisateur, il est dans le répertoire `/data/omics-school/share/miniconda/envs/rnaseq`
et est activable avec la commande
```
$ conda activate rnaseq
```

Un autre utilisateur (par exemple `adejardin`) ne peut pas modifier l'environnement `rnaseq` ni même créer un nouvel environnement dans le répertoire par défaut de conda (`/data/omics-school/share/miniconda/envs/`) car il n'en pas les droits.

La solution est alors de créer un nouvel environnement en indiquant dans quel répertoire cet environnement doit être créé. Par exemple :
```
$ conda create -y -p $HOME/conda-env-projet
```

La commande pour activer l'environnement est :
```
$ conda activate $HOME/conda-env-projet
```

Pour l'utilisateur `adejardin`, l'invite de commande doit alors ressembler à quelque chose du type :
```
(/data/omics-school/adejardin/env-conda-env-projet) adejardin@candihub:~$
```

Il faut ensuite installer les logiciels nécessaires, par exemple :
```
$ conda install -y fastqc star htseq samtools
```

Rappel : pour quitter l'environnement une fois que les analyses sont terminées :
```
$ conda deactivate
```
