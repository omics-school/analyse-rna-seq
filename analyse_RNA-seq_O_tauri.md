---
title: Analyse des données RNA-seq de O. tauri sous Linux
author: Pierre Poulain
license: Creative Commons Attribution (CC-BY)
---

# Analyse des données RNA-seq de O. tauri sous Linux

Dans cette activité, nous allons analyser les données RNA-seq de O. tauri sous Linux.

Pour cela, vous allez beaucoup utiliser la ligne de commande, connecté en SSH sur le serveur du DU. Vous copierez également des fichiers depuis le serveur du DU vers votre machine locale (avec FileZilla).

Voici une vue d'ensemble des étapes pour analyser les données de séquençage haut débit de *O. tauri.* :

![](pipeline_RNA_seq_O_tauri.svg)


## Préparation de l'environnement de travail

Pour analyser les données de séquençage haut débit de *O. tauri.*, nous avons besoin des logiciels suivants :

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) pour le contrôle qualité des données de séquençage.
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) pour l'indexation du génome de référence puis l'alignement des *reads* sur le génome de référence.
- [SAMtools](http://samtools.sourceforge.net/) pour la manipulation des fichiers d'alignements (conversion en binaire, tri et indexation)
- [HTSeq](https://htseq.readthedocs.io/en/latest/) pour le comptage du nombre de *reads* alignés sur chaque gène.


### Anaconda, Miniconda, Conda & Bioconda

[Anaconda](https://www.anaconda.com/what-is-anaconda/) est une distribution open source, disponible pour Windows, Mac et Linux, et qui contient de nombreux outils utilisés pour l'analyse de données avec le langage de programmation Python.

Anaconda est basé sur [conda](https://conda.io/docs/), un gestionnaire de paquets, qui permet d'installer des logiciels facilement et sans être administrateur. Conda permet d'installer des logiciels écrits en Python mais aussi en R, en C...

Enfin, Anaconda est également disponible dans une version *light* appelée [Miniconda](https://conda.io/miniconda.html). Miniconda ne contient pas tous les outils Python disponibles dans Anaconda, mais il contient néanmoins le gestionnaire de paquets conda.

Enfin, [Bioconda](https://bioconda.github.io/) est un canal de diffusion de logiciels, utilisable par le gestionnaire de paquets conda et proposant de nombreux logiciels utilisés en bioinformatique.

Voici deux articles très intéressants sur conda :
- [Conda le meilleur ami du bioinformaticien](https://bioinfo-fr.net/conda-le-meilleur-ami-du-bioinformaticien). Article d'introduction. Attention cependant, certaines commandes sont obsolètes.
- [Comment fixer les problèmes de déploiement et de durabilité des outils en bioinformatique ? Indice : conda !](https://bioinfo-fr.net/comment-fixer-les-problemes-de-deploiement-et-de-durabilite-des-outils-en-bioinformatique). Article un peu plus technique.


### Configuration de conda

La distribution Miniconda a été installée sur le serveur du DU. Pour que vous puissiez avoir accès à cet outil,
vous devez configurer votre *shell* Linux sur le serveur du DU. Les étapes à suivre sont :

1. Connectez-vous en SSH au serveur du DU.
1. Éditez le fichier `.bashrc` dans votre répertoire personnel. Par exemple avec l'éditeur de texte nano :
    ```
    $ nano .bashrc
    ```
1. Déplacez-vous à la fin du fichier et ajoutez la ligne ci-dessous :
    ```
    source /data/omics-school/share/miniconda/etc/profile.d/conda.sh
    ```
    Enregistrez le fichier (`ctrl + o`) puis quittez nano (`ctrl + x`).  
    Remarque 1 : la ligne de commande à ajouter est assez longue. Pour éviter les erreurs, utilisez le copier / coller (`ctrl + maj + v`) dans nano.  
    Remarque 2 : il est possible que votre fichier `.bashrc` soit vide, ce n'est pas un problème.
1. Vérifiez que conda est maintenant disponible en tapant les commandes suivantes :
    ```
    $ source .bashrc
    $ conda --version
    ```
1. Bravo ! :tada: Vous avez correctement configuré conda. Fermez votre shell sur le serveur du DU.

Les manipulations ci-dessus vous ont permis de rendre disponible conda dans votre *shell* Linux sur le serveur du DU. Elles ne sont à faire qu'une seule fois.

Une documentation expliquant l'installation de miniconda et la configuration de conda sur le serveur du DU est disponible [ici](installation_conda_logiciels_RNA-seq.md).


### Chargement de l'environnement conda

Conda est un gestionnaire de paquets qui permet d'installer de nombreux logiciels utilisés en bioinformatique. Il permet aussi de créer des environnements *virtuels* dans lequel ces logiciels sont installés. L'intérêt des environnements virtuels est de pouvoir installer sur la même machine, plusieurs version d'un même logiciel, chaque version étant installée dans un environnement virtuel indépendant.

Nous allons maintenant voir comment charger un environnement virtuel créé avec conda.

1. Connectez-vous en SSH au serveur du DU.
1. Vérifiez que conda est bien disponible avec la commande
    ```
    $ conda --version
    ```
1. Nous avons préparé un environnement virtuel conda spécialement pour l'analyse RNA-seq. Chargez cet environnement :
    ```
    $ conda activate rnaseq
    ```
1. Votre invite de commande devrait être modifiée et ressembler à :
    ```
    (rnaseq) ppoulain@candihub:~$
    ```
    La mention `(rnaseq)` indique que vous êtes dans l'environnement virtuel `rnaseq`.

Remarque : pour quitter un environnement virtuel, il faut utiliser la commande
```
$ conda deactivate
```
mais vous n'aurez pas besoin de l'utiliser pour cette activité.

Pour la suite, nous supposerons que :
1. Vous êtes connecté en SSH au serveur du DU.
1. Vous avez activé l'environnement conda `rnaseq`.

### Vérification des logiciels disponibles

En bioinformatique, il est essentiel de vérifier et noter les versions des logiciels que vous utilisez.

Dans votre environnement virtuel conda, vérifiez les versions des logiciels que vous allez utiliser :

```
$ fastqc --version
$ bowtie2 --version | head -n 1
$ samtools --version | head -n 1
$ htseq-count -h | grep version
```

Si vous vous demandez pourquoi on utilise parfois `| head -n 1` ou `| grep version`, comparez par exemple les sorties des commandes
```
$ fastqc --version
```
et
```
$ bowtie2 --version
```
Certains programmes peuvent renvoyer beaucoup d'informations.


### Comparaison avec les logiciels utilisés dans Galaxy

Connectez-vous maintenant à votre compte sur Galaxy. Essayez de retrouver les versions des logiciels que vous utilisés (FastQC, Bowtie2, SAMtools, HTSeq).

Pour ce faire, dans votre *History*, cliquez sur le nom d'un résultat d'analyse, puis cliquez sur le petit i entouré (:information_source:) et lisez les informations de la section *Job Dependencies*.


## Préparation des données

Sur le serveur du DU, dans votre répertoire personnel, créez le répertoire `RNAseq`.

Dans ce répertoire, copiez :

- Les 2 ou 3 fichiers contenant les *reads* (`.fastq.gz`) qui vous ont été attribués la dernière fois et que vous avez analysés avec Galaxy. Tous les fichiers sont dans le répertoire  `/data/omics-school/share/RNAseq_tauri/`
- Le génome de référence de *O. tauri* : `/data/omics-school/share/GCF_000214015.3_version_140606_genomic.fna`
- Les annotations du génome de référence :
`/data/omics-school/share/GCF_000214015.3_version_140606_genomic_DUO2.gff`

Remarque : le génome de référence de *Ostreococcus tauri* et ses annotations sont disponibles sur la [page dédiée](https://www.ncbi.nlm.nih.gov/genome/373?genome_assembly_id=352933) sur le site du NCBI :
- [génome de référence](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.fna.gz)
- [annotations](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.gff.gz). Nous avons légèrement modifié le fichier d'annotations pour ne prendre en compte que les gènes et alléger la visualisation dans IGV.


**:warning: Étape essentielle pour la suite :warning:**. Renommez les fichiers contenant vos *reads* (`.fastq.gz`) de la façon suivante :
```
HCA-XX_R1.fastq.gz
```
avec `XX` le numéro de votre échantillon. Cela revient à supprimer la première partie du nom du fichier.

Par exemple, pour les fichiers correspondants aux échantillons 10 et 41 :
```
$ mv 140317_SN365_A_L001_HCA-10_R1.fastq.gz HCA-10_R1.fastq.gz
$ mv 140317_SN365_A_L002_HCA-41_R1.fastq.gz HCA-41_R1.fastq.gz
```

Cela simplifiera considérablement l'automatisation des analyses par la suite.

## Analyse manuelle

Pour cette première analyse, choisissez un **seul échantillon** contenant des *reads*.


### Contrôle qualité

Lancez FastQC avec la commande :

```
$ fastqc nom-fichier-fastq.gz
```
où `nom-fichier-fastq.gz` est le fichier contenant l'échantillon que vous avez choisi.

FastQC va produire deux fichiers (`.html` et `.zip`). Copiez le fichier `.html` sur votre machine locale avec FileZilla.

Ouvrez ce fichier dans un navigateur internet (Firefox par exemple).

Analysez le rapport de FastQC.


### Indexation du génome de référence

Sur le serveur du DU, dans le répertoire `RNAseq`, lancez l'indexation du génome de référence.
```
$ bowtie2-build GCF_000214015.3_version_140606_genomic.fna O_tauri
```
Les index vont être stockés dans des fichiers dont le nom débute par `O_tauri`.

Calculez la taille total des fichiers index avec la commande
```
$ du -ch O_tauri*
```

Explications : la commande `du` affiche la taille occupée par des fichiers. L'option `-h` affiche la taille en ko, Mo, Go... L'option `-c` calcule la taille totale occupée par tous les fichiers.


### Alignements des *reads* sur le génome de référence

Lancez l'alignement :
```
$ bowtie2 -x O_tauri -U nom-fichier-fastq.gz -S bowtie.sam
```

Ici :
- `O_tauri` désigne les fichiers index du génome de référence,
- `nom-fichier-fastq.gz` est le fichier contenant l'échantillon que vous avez choisi
- et `bowtie.sam` est le fichier qui va contenir l'alignement produit par Bowtie2.

Cette étape est la plus longue et peut prendre plusieurs minutes (~ 10).

À la fin de l'alignement, Bowtie2 renvoie plusieurs informations intéressantes comme :
- le nombre total de *reads* lus dans le fichier `.fastq.gz`
- le nombre de *reads* non alignés (*aligned 0 times*)
- le nombre de *reads* alignés une seule fois
- le nombre de *reads* alignés plus d'une fois
- un taux d'alignement global

Il faut être prudent si le nombre de *reads* non alignés est trop important (> 20%).


### Conversion des *reads* alignés en binaire, tri et indexation

Vous allez maintenant utiliser SAMtools pour :

1. Convertir le fichier `.sam` créé par Bowtie2, qui est un fichier texte, en fichier `.bam`, qui est un fichier binaire, et qui donc prend moins de place sur le disque.
    ```
    $ samtools view -b bowtie.sam > bowtie.bam
    ```
1. Trier les *reads* alignés suivant l'ordre dans lequel ils apparaissent dans le génome.
    ```
    $ samtools sort bowtie.bam -o bowtie.sorted.bam
    ```
1. Indexer le fichier `.bam`. Cette étape est indispensable pour visualiser l'alignement avec IGV.
    ```
    $ samtools index bowtie.sorted.bam
    ```


### Visualisation des *reads* alignés avec IGV

Pour visualiser l'alignement des *reads* sur le génome de référence avec IGV, copiez, avec FileZilla, sur votre machine locale les fichiers :
- génome de référence (fichier `.fna`) ;
- annotations du génome de référence (fichier `_DUO2.gff`) ;
- bam trié (`bowtie.sorted.bam`) ;
- index du bam trié (`bowtie.sorted.bam.bai`).

Lancez IGV et visualisez l'alignement des *reads* sur le génome de référence. Si vous avez oublié comme faire, visionnez la vidéo 2, de l'activité 1, du cours de Mai sur cloudschool.

Visualisez particulièrement le gène `ostta18g01980`.


### Comptage des *reads* alignés sur les gènes de *O. tauri*

Le comptage des *reads* alignés sur les gènes se fait avec HTSeq.

Lancez la commande :
```
$ htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam bowtie.sorted.bam GCF_000214015.3_version_140606_genomic_DUO2.gff > count.txt
```

HTSeq renvoie le nombre d'annotations trouvées dans le fichier `.gff`.

Déterminez le nombre de *reads* alignés sur le gène `ostta18g01980`. Pour cela, vous pouvez lancer la commande
```
$ grep ostta18g01980 count.txt
```
ou alors ouvrir le fichier `count.txt` avec la commande `less` puis chercher `ostta18g01980` en tapant `/ostta18g01980` puis la touche `Entrée`.


## Automatisation de l'analyse : niveau 1

Tout cela est très bien mais les fichiers que vous avez générés (`bowtie.bam`, `bowtie.sorted.bam`, `count.txt`...) ne sont pas très informatifs sur l'échantillon dont ils proviennent.

Par ailleurs, entrer toutes ces commandes à la main, les unes après les autres, est pénible et source d'erreurs.

Pour répondre à ces deux problèmes, nous allons introduire les notions Bash de variables et de  scripts.

### Variables

Une variable va simplement contenir de l'information qui sera utilisable autant de fois que nécessaire.

Création de variables :
```
$ toto=33
$ t="salut"
```
Il faut coller le nom de la variable et son contenu au symbole `=`.

Affichage de variables :
```
$ echo $toto
33
$ echo "$t Pierre"
salut Pierre
```
La commande `echo` permet d'afficher une chaîne de caractère, une variable, ou les deux.

Pour utiliser une variable (et accéder à son contenu), il faut précéder son nom du caractère `$`. Attention, ce symbole n'est pas à confondre avec celui qui désigne l'invite de commande de votre *shell* Linux.

Enfin, une bonne pratique consiste à utiliser une variable avec le symbole `$` et son nom entre accolades :
```
$ echo ${toto}
33
$ echo "${t} Pierre"
salut Pierre
```

### Script

Un script est un fichier texte qui contient des instructions Bash. Par convention, il porte l'extension `.sh`.

Dans un script Bash, tout ce qui suit le symbole `#` est considéré comme un commentaire et n'est donc pas traité par Bash.


### Analyse RNA-seq

Observez le script bash [script1.sh](script1.sh) et essayer de comprendre son fonctionnement, notamment l'utilisation des variables.

Testez le script `script1.sh` sur **un seul** de vos échantillons. Pour cela :
- Recopiez le script dans un fichier `script1.sh` dans votre répertoire `RNAseq` ou, plus simplement, téléchargez-le directement avec la commande
```
$ wget https://raw.githubusercontent.com/omics-school/analyses-rna-seq-o-tauri/master/script1.sh
```
- Ouvrez le script `script1.sh` avec `nano` et modifiez la variable `sample` avec votre numéro d'échantillon. Sauvegardez le script (`ctrl + o`) et quittez nano (`ctrl + x`).  
Rappel : pas d'espace avant ou après le symbole `=` !
- Lancez le script avec la commande
    ```
    $ bash script1.sh
    ```

Vérifiez que le déroulement du script se passe bien. Vous avez le temps de prendre un café :coffee:. Voir plusieurs :coffee: :cookie: :coffee: :cookie:.


## Automatisation de l'analyse : niveau 2

Le script précédent était pratique mais il ne conserve pas les informations liées à l'alignement (nombre de *reads* non-alignés, alignés une fois...).

Proposez une évolution du premier script pour répondre à ce problème. N'hésitez pas à améliorer aussi la lisibilité du script lors de son exécution.

La solution est dans le script 2. Pour le télécharger, utilisez la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyses-rna-seq-o-tauri/master/script2.sh
```
Vous remarquerez que la solution proposée pour conserver les informations liées à l'alignement est un peu particulière. Nous allons en discuter, mais dans un premier temps essayer de comprendre l'explication donnée [ici](https://stackoverflow.com/questions/876239/how-can-i-redirect-and-append-both-stdout-and-stderr-to-a-file-with-bash).


## Automatisation de l'analyse : niveau 3 (ninja)

Le script précédent était intéressant mais il ne prend en compte qu'un seul échantillon à la fois.

On aimerait avoir un seul script qui traiterait tous les échantillons qu'on souhaite analyser.
Cela est possible avec une boucle. Une boucle permet de répéter un ensemble d'instructions.

Voici un exemple en Bash :
```
$ for prenom in gaelle bertrand pierre
> do
> echo "Salut ${prenom} !"
> done
Salut gaelle !
Salut bertrand !
Salut pierre !
```
En sacrifiant un peu de lisibilité, la même commande peut s'écrire sur une seule ligne :
```
$ for prenom in gaelle bertrand pierre; do echo "Salut ${prenom} !"; done
Salut gaelle !
Salut bertrand !
Salut pierre !
```

Notez l'utilisation du symbole `;` pour séparer les différents éléments de la boucle.

Une leçon de Software Carpentry aborde la notion de [boucle](https://swcarpentry.github.io/shell-novice/05-loop/index.html). Prenez quelques minutes pour la parcourir et faire les exercices.

Le [script 3](script3.sh) utilise une boucle. Observez la structure du script et essayez de comprendre son fonctionnement.

La ligne `set -euo pipefail` en tout début du script va arrêter celui-ci :
- à la première erreur
- si une variable n'est pas définie
- si une erreur est rencontrée dans une commande avec un pipe (`|`)

C'est une mesure de sécurité importante pour votre script. Si vous le souhaitez, vous pouvez lire l'article de Aaron Maxwell à ce sujet : [Use the Unofficial Bash Strict Mode (Unless You Looove Debugging)](http://redsymbol.net/articles/unofficial-bash-strict-mode/)

Téléchargez le script 3 avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyses-rna-seq-o-tauri/master/script3.sh
```

Modifiez le script 3 avec les numéros d'échantillons que vous avez à analyser. Faites bien attention à la variable concernée et sa syntaxe.

Si vous pensez en avoir le temps, lancez le script 3.

**Remarque**. Au tout début de l'activité, vous avez renommé les fichiers contenant les *reads* (`.fastq.gz`). Cette étape manuelle peut être automatisée, par exemple avec la commande suivante :

```
for name in $(ls 140317*.fastq.gz); do newname=$(echo ${name} | sed 's/140317_SN365_A_L00._//'); mv ${name} ${newname}; done
```

Essayez de comprendre son fonctionnement.
