---
title: Analyse des données RNA-seq
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

L'objectif de cette partie est d'analyser les données RNA-seq de *O. tauri* sur votre machine Unix locale.

Voici une vue d'ensemble des étapes pour analyser ces données de séquençage haut débit :

![](pipeline_RNA_seq_O_tauri.svg)


## 3.1 Préparer l'environnement de travail

Sous Windows, ouvrez un terminal Ubuntu.

Déplacez-vous dans le répertoire `/mnt/c/Users/omics/rnaseq_tauri` :
```
$ cd /mnt/c/Users/omics/rnaseq_tauri
```

Activez l'environnement conda *rnaseq-env* :
```
$ conda activate rnaseq-env
```

*Remarque : contrôlez que le nom de l'environnement conda apparait bien à gauche de l'invite de commande sous la forme : `(rnaseq-env)`.*

Vous êtes maintenant prêt à analyser des données RNA-seq 🤠


## 3.2 Analyse manuelle

Pour cette première analyse, choisissez un **seul échantillon** contenant des *reads*, c'est-à-dire un fichier parmi :
```
reads/SRR2960338.fastq.gz
reads/SRR2960341.fastq.gz
reads/SRR2960343.fastq.gz
```

### 3.2.1 Contrôler la qualité des reads

Créez le répertoire `reads_qc` qui va contenir les fichiers produits par le contrôle qualité des fichiers fastq.gz :

```bash
$ mkdir -p reads_qc
```

Lancez FastQC avec la commande :

```bash
$ fastqc reads/ECHANTILLON.fastq.gz --outdir reads_qc
```
où `reads/ECHANTILLON.fastq.gz` est le fichier contenant l'échantillon que vous avez choisi. Pensez à l'adapter avant de lancer votre commande !

FastQC va produire deux fichiers (un fichier avec l'extension `.html` et un autre avec l'extension `.zip`) dans le répertoire `reads_qc`. Si par exemple, vous avez analysé le fichier `reads/SRR2960338.fastq.gz`, vous obtiendrez les fichiers `reads_qc/SRR2960338_fastqc.html` et `reads_qc/SRR2960338_fastqc.zip`.

Utilisez l'explorateur de fichiers de Windows pour vous déplacer dans le bon répertoire (voir [exemple](img/navigation_explorateur_windows.png)). Ouvrez ensuite le fichier `.html` ainsi créé avec Firefox (en double-cliquant sur le fichier). Analysez le rapport créé par FastQC.


### 3.2.2 Indexer le génome de référence

L'indexation du génome de référence est une étape indispensable pour accélérer l'alignement des reads sur le génome. Elle consiste à créer un annuaire du génome de référence.

Toujours depuis votre shell Ubuntu et dans le répertoire `/mnt/c/Users/omics/rnaseq_tauri`, créez le répertoire `index` :

```bash
$ mkdir -p index
```

Lancez l'indexation du génome de référence.

```bash
$ bowtie2-build genome/GCF_000214015.3_version_140606.fna index/O_tauri
```

Les index sont stockés dans des fichiers dont le nom débute par `O_tauri` dans le répertoire `index`.

Calculez la taille total des fichiers index avec la commande :

```bash
$ du -ch index/O_tauri*
```

Comparez la taille totale des index à la taille du fichier contenant le genome (`genome/GCF_000214015.3_version_140606.fna`).

L'indexation du génome n'est à faire qu'une seule fois pour chaque logiciel d'alignement.


### 3.2.3 Aligner les *reads* sur le génome de référence

Créez le répertoire `map` qui va contenir les *reads* alignés sur le génome de référence :

```bash
$ mkdir -p map
```

Lancez l'alignement :

```bash
$ bowtie2 -p 2 -x index/O_tauri -U reads/nom-fichier.fastq.gz -S map/bowtie.sam
```

Les options utilisées sont :

- `-p 2` : utilisation de 2 coeurs pour réaliser l'alignement. Votre machine possède 4 coeurs.
- `-x index/O_tauri` : désigne les fichiers index du génome de référence.
- `-U reads/ECHANTILLON.fastq.gz` : indique le nom du fichier contenant les reads. Adaptez-le au nom de l'échantillon que vous avez choisi.
- `-S map/bowtie.sam` : précise le nom du fichier de sortie qui va contenir l'alignement produit par Bowtie2.

Toutes les options de Bowtie2 sont détaillées dans la [documentation](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml).

Cette étape peut prendre 6 ou 7 minutes. **Bowtie n'affiche rien à l'écran lorsqu'il fonctionne**. Soyez patient.

À la fin de l'alignement, Bowtie2 renvoie des informations qui ressemblent à :

```
5470272 reads; of these:
  5470272 (100.00%) were unpaired; of these:
    497809 (9.10%) aligned 0 times
    4603112 (84.15%) aligned exactly 1 time
    369351 (6.75%) aligned >1 times
90.90% overall alignment rate
```

On obtient ainsi :
- le nombre total de *reads* lus dans le fichier `.fastq.gz` (ici `5470272`)
- le nombre de *reads* non alignés « *aligned 0 times* » (`497809`, soit `9.10%` du nombre total de *reads*)
- le nombre de *reads* alignés une seule fois (`4603112`, soit `84.15%` du nombre total de *reads*)
- le nombre de *reads* alignés plus d'une fois (`369351`, soit `6.75%` du nombre total de *reads*)
- un taux d'alignement global (`90.90%`)

Il faut être prudent si le nombre de *reads* non alignés est trop important (par exemple supérieur à 20%).


### 3.2.4 Convertir en binaire, trier et indexer les *reads* alignés

Vous allez maintenant utiliser SAMtools pour :

1. Convertir le fichier `.sam` (fichier texte) créé par Bowtie2 en fichier `.bam`, qui est un fichier binaire, et qui donc prend moins de place sur le disque.
    ```bash
    $ samtools view -@ 2 -b map/bowtie.sam -o map/bowtie.bam
    ```
    Cette étape va prendre plusieurs minutes. Comme votre machine dispose de 4 coeurs, nous allons en utiliser 2 (`-@ 2`) pour accélérer le calcul.  
    Avec la commande `du` et les bonnes options, comparez la taille deux fichiers `map/bowtie.sam` et `map/bowtie.bam`. Quel est le ratio de compression entre les deux formats de fichiers ?

2. Trier les *reads* alignés suivant l'ordre dans lequel ils apparaissent dans le génome.
    ```bash
    $ samtools sort -@ 2 map/bowtie.bam -o map/bowtie.sorted.bam
    ```
    Cette étape va prendre ici encore quelques minutes.

3. Indexer le fichier `.bam`. Cette étape est **indispensable** pour visualiser l'alignement avec IGV.
    ```bash
    $ samtools index map/bowtie.sorted.bam
    ```


### 3.2.5 Compter les *reads* alignés sur les gènes de *O. tauri*

Le comptage des *reads* alignés sur les gènes se fait avec HTSeq.

Toujours depuis votre shell Ubuntu et dans le répertoire `/mnt/c/Users/omics/rnaseq_tauri`, créez le répertoire `count` :

```bash
$ mkdir -p count
```

Puis lancez la commande (en une seule ligne) pour compter les *reads* alignés :

```bash
$ htseq-count --stranded=no --type="gene" --idattr="ID" --order=name --format=bam map/bowtie.sorted.bam genome/GCF_000214015.3_version_140606.gff > count/count.txt
```

HTSeq renvoie le nombre d'annotations trouvées dans le fichier `.gff` (7659) puis affiche une progression de l'analyse. Les options du programme `htseq-count` sont décrites dans la [documentation](https://htseq.readthedocs.io/en/master/htseqcount.html).

Déterminez le nombre de *reads* alignés sur le gène `ostta18g01980` avec la commande :

```bash
$ grep ostta18g01980 count/count.txt
```

Vous pouvez également ouvrir le fichier `count/count.txt` avec la commande `less` puis chercher le gène `ostta18g01980` en tapant `/ostta18g01980` puis la touche <kbd>Entrée</kbd> (et enfin la touche <kbd>Q</kbd> pour quitter).


### 3.2.6 Visualiser les *reads* alignés avec IGV

Pour visualiser l'alignement des *reads* sur le génome de référence avec IGV, vous avez besoin des fichiers suivants :
- Le génome de référence (`genome/GCF_000214015.3_version_140606.fna`).
- Les annotations du génome de référence (`genome/GCF_000214015.3_version_140606.gff`).
- Le fichier bam trié (`map/bowtie.sorted.bam`).
- L'index du bam trié (`map/bowtie.sorted.bam.bai`).

Lancez IGV et visualisez l'alignement des *reads* sur le génome de référence. Si vous avez oublié comme faire, visionnez la vidéo sur ce sujet qui vous a été proposée précédemment.

Visualisez particulièrement le gène `ostta18g01980`.


## 3.3 Automatiser l'analyse : niveau 1

Tout cela est très bien mais les fichiers que vous avez générés (`map/bowtie.bam`, `map/bowtie.sorted.bam`, `cout/count.txt`...) portent des noms qui ne sont pas très informatifs sur l'échantillon dont ils proviennent.

Par ailleurs, entrer toutes ces commandes à la main, les unes après les autres, est pénible et source d'erreurs. Et il y a fort à parier que vous aurez complètement oublié ces commandes dans 1 semaine, voire dans 1 heure. 🤯 C'est parfaitement normal, il n'y a absolument aucun intérêt à se souvenir de toutes ces commandes.

Pour répondre à ces deux problèmes, de gestion de données et d'automatisation, nous allons introduire les notions Bash de variables et de scripts.

Mais d'abord, faites un peu de ménage en supprimant les fichiers créés précédemment :

```bash
$ rm -f reads_qc/*fastqc* index/*bt2 map/bowtie* count/count*
```

*💣 Attention à l'utilisation de la commande `rm` qui supprime définitivement les fichiers. 💣*


### 3.3.1 Variables

Une variable va simplement contenir de l'information qui sera utilisable autant de fois que nécessaire.

Création de variables :

```bash
$ toto=33
$ t="salut"
```

*Attention : Il faut coller le nom de la variable et son contenu au symbole `=`.*

Affichage de variables :

```bash
$ echo $toto
33
$ echo "$t Pierre"
salut Pierre
```

La commande `echo` affiche une chaîne de caractère, une variable, ou les deux.

Pour utiliser une variable (et accéder à son contenu), il faut précéder son nom du caractère `$`. Attention, ce symbole n'est pas à confondre avec celui qui désigne l'invite de commande de votre *shell* Linux.

Enfin, une bonne pratique consiste à utiliser une variable avec le symbole `$` et son nom entre accolades :

```bash
$ echo ${toto}
33
$ echo "${t} Pierre"
salut Pierre
```

### 3.3.2 Script

Un script est un fichier texte qui contient des instructions Bash. Par convention, il porte l'extension `.sh`. L'objectif premier d'un script Bash est d'automatisr l'exécution de plusieurs commandes Bash, la plupart du temps pour manipuler ou analyser des fichiers.

Dans un script Bash, tout ce qui suit le symbole `#` est considéré comme un commentaire et n'est donc pas traité par Bash.


### 3.3.3 Analyse RNA-seq

Testez le script `script1.sh` sur **un seul** de vos échantillons. Pour cela :

- Vérifiez que vous êtes bien toujours dans le répertoire `/mnt/c/Users/omics/rnaseq_tauri`.

- Téléchargez le script `script1.sh` dans votre répertoire `rnaseq_tauri` avec la commande :
    ```bash
    $ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script1.sh
    ```

- Ouvrez le script `script1.sh` avec `nano`. Essayez de comprendre son fonctionnement, notamment l'utilisation des variables.

- Sur la deuxième ligne, modifiez la variable `sample` avec votre numéro d'échantillon. Sauvegardez le script (<kbd>Ctrl</kbd> + <kbd>O</kbd>) et quittez nano (<kbd>Ctrl</kbd> + <kbd>X</kbd>).  
    Rappel : pas d'espace avant et après le symbole `=` !

- Lancez le script avec la commande :
    ```bash
    $ bash script1.sh
    ```

Vérifiez que le déroulement du script se passe bien. Vous avez le temps de prendre un café (~ 15 '), voir plusieurs ☕ 🍪 ☕ 🍪.

Évaluez approximativement le temps nécessaire au script 1 pour s'exécuter. ⏱️ À partir de cette valeur, extrapoler le temps nécessaire qu'il faudrait pour analyser les 3 échantillons.

Utilisez enfin la commande `tree` pour contempler votre travail (ici avec l'échantillon SRR2960338) :

```bash
$ tree
.
├── count
│   └── count-SRR2960338.txt
├── genome
│   ├── GCF_000214015.3_version_140606.fna
│   ├── GCF_000214015.3_version_140606.gff
│   └── md5sum.txt
├── index
│   ├── O_tauri.1.bt2
│   ├── O_tauri.2.bt2
│   ├── O_tauri.3.bt2
│   ├── O_tauri.4.bt2
│   ├── O_tauri.rev.1.bt2
│   └── O_tauri.rev.2.bt2
├── map
│   ├── bowtie-SRR2960338.sorted.bam
│   └── bowtie-SRR2960338.sorted.bam.bai
├── reads
│   ├── SRR2960338.fastq.gz
│   ├── SRR2960341.fastq.gz
│   └── SRR2960343.fastq.gz
├── reads_qc
│   ├── SRR2960338_fastqc.html
│   └── SRR2960338_fastqc.zip
└── script1.sh

6 directories, 18 files
```


## 3.4 Automatiser l'analyse : niveau 2

Le script précédent est pratique mais il ne conserve pas les informations liées à l'alignement générées par Bowtie2 (nombre de *reads* non-alignés, alignés une fois...).

Le [script 2](script2.sh) répond à ce problème. Téléchargez-le avec la commande :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script2.sh
```

Ouvrez ce script avec `nano`. Vous remarquerez que la solution proposée pour conserver les informations liées à l'alignement est un peu particulière (`2> map/bowtie-${sample}.out`). Nous allons en discuter, mais dans un premier temps essayer de comprendre l'explication donnée [ici](https://stackoverflow.com/questions/876239/how-can-i-redirect-and-append-both-stdout-and-stderr-to-a-file-with-bash).


## 3.5 Automatiser l'analyse : niveau 3 (ninja)

Le script précédent est intéressant mais il ne prend en compte qu'un seul échantillon à la fois. Quel ennui !

On aimerait avoir un seul script qui traiterait tous les échantillons qu'on souhaite analyser.
Cela est possible avec une boucle. Une boucle permet de répéter un ensemble d'instructions.

Voici un exemple en Bash :

```bash
$ for prenom in gaelle bertrand pierre
> do
> echo "Salut ${prenom} !"
> done
Salut gaelle !
Salut bertrand !
Salut pierre !
```

En sacrifiant un peu de lisibilité, la même commande peut s'écrire sur une seule ligne :

```bash
$ for prenom in gaelle bertrand pierre; do echo "Salut ${prenom} !"; done
Salut gaelle !
Salut bertrand !
Salut pierre !
```

Notez l'utilisation du symbole `;` pour séparer les différents éléments de la boucle.

Une leçon de Software Carpentry aborde la notion de [boucle](https://swcarpentry.github.io/shell-novice/05-loop/index.html). Prenez quelques minutes pour parcourir cette leçon et comprendre de quoi il s'agit.

Le [script 3](script3.sh) utilise une boucle pour automatiser l'analyse de plusieurs échantillons. Téléchargez-le avec la commande :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script3.sh
```

Ouvrez ce script avec `nano`. Observez la structure du script et essayez de comprendre son fonctionnement.

La ligne `set -euo pipefail` tout au début du script va arrêter celui-ci :
- à la première erreur ;
- si une variable n'est pas définie ;
- si une erreur est rencontrée dans une commande avec un pipe (`|`).

C'est une mesure de sécurité importante pour votre script. Si vous le souhaitez, vous pouvez lire l'article de Aaron Maxwell à ce sujet : [Use the Unofficial Bash Strict Mode (Unless You Looove Debugging)](http://redsymbol.net/articles/unofficial-bash-strict-mode/)

Si vous pensez en avoir le temps, lancez le script 3. Comme ce script va automatiser toute l'analyse, il va fonctionner environ 45 minutes et vous aurez peut-être besoin de fermer votre terminal. Pour ne pas arrêter brutalement l'analyse, lancez le script de cette manière :

```bash
$ nohup bash script3.sh &
```

Mais pour autant, n'éteignez pas votre ordinateur !

Le message 

```bash
nohup: ignoring input and appending output to 'nohup.out'
```
vous rappelle que les messages qui apparaissaient habituellement à l'écran seront redirigés dans le fichier `nohup.out`.


##  3.6 Comparer les versions des logiciels utilisés dans Galaxy (si vous avez du temps)

Connectez-vous maintenant à votre compte sur Galaxy. Essayez de retrouver les versions des logiciels que vous avez utilisés (FastQC, Bowtie2, SAMtools, HTSeq).

Pour ce faire, dans votre *History*, cliquez sur le nom d'un résultat d'analyse, puis cliquez sur le petit i entouré (ℹ️) et lisez les informations de la section *Job Dependencies*.

Comparez les versions des logiciels disponibles dans Galaxy et de ceux que vous avez utilisés sur votre machine.
