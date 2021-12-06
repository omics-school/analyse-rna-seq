---
title: Analyse des données RNA-seq
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

L'objectif de cette partie est d'analyser les données RNA-seq de *O. tauri* sur votre machine Unix locale.

Voici une vue d'ensemble des étapes pour analyser les données de séquençage haut débit :

![](pipeline_RNA_seq_O_tauri.svg)


## 3.1 Préparation de l'environnement de travail

Sous Windows, ouvrez un terminal Ubuntu.

Déplacez-vous dans le répertoire `/mnt/c/Users/omics/rnaseq_tauri` :
```
$ cd /mnt/c/Users/omics/rnaseq_tauri
```

Activez l'environnement conda *rnaseq-env* :
```
$ conda activate rnaseq-env
```

Remarque : contrôlez que le nom de l'environnement conda apparait bien à gauche de l'invite de commande sous la forme : `(rnaseq-env)`.

Vous êtes maintenant prêt à analyser des données RNA-seq 🤠


## 3.2 Analyse manuelle

Pour cette première analyse, choisissez un **seul échantillon** contenant des *reads*, c'est-à-dire un fichier parmi :
```
reads/SRR2960338.fastq.gz
reads/SRR2960341.fastq.gz
reads/SRR2960343.fastq.gz
```

### Contrôle qualité

Créez le répertoire `reads_qc` qui va contenir les fichiers produits par le contrôle qualité des fichiers fastq.gz :

```bash
$ mkdir -f reads_qc
```

Lancez FastQC avec la commande :

```bash
$ fastqc reads/nom-fichier.fastq.gz --outdir reads_qc
```
où `nom-fichier.fastq.gz` est le fichier contenant l'échantillon que vous avez choisi.

FastQC va produire deux fichiers (un fichier avec l'extension `.html` et un autre avec l'extension `.zip`) dans le répertoire `reads_qc`. Si par exemple, vous avez analysé le fichier `reads/SRR2960338.fastq.gz`, vous obtiendrez les fichiers `reads_qc/SRR2960338_fastqc.html` et `reads_qc/SRR2960338_fastqc.zip`.

En utilisant l'explorateur de fichiers de Windows, ouvrez le fichier `.html` ainsi créé avec Firefox (en cliquant sur le fichier). Analysez le rapport de FastQC.


### Indexation du génome de référence

L'indexation du génome de référence est une étape indispensable pour accélérer l'alignement des reads sur le génome.

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


### Alignements des *reads* sur le génome de référence

Créez le répertoire `map` qui va contenir les *reads* alignés sur le génome de référence :
```
$ mkdir -p map
```

Lancez l'alignement :
```
$ bowtie2 -p 2 -x index/O_tauri -U reads/nom-fichier.fastq.gz -S map/bowtie.sam
```

Ici :
- `genome/O_tauri` désigne les fichiers index du génome de référence,
- `reads/nom-fichier.fastq.gz` est le fichier contenant l'échantillon. Adaptez-le au nom de l'échantillon que vous avez choisi.
- et `reads/bowtie.sam` est le fichier qui va contenir l'alignement produit par Bowtie2.

Comme votre machine dispose de 4 coeurs, nous en utilisons 2 (`-p 2`) pour accélérer le calcul.

Cette étape peut prendre plusieurs minutes. **Bowtie n'affiche rien à l'écran lorsqu'il fonctionne**. Soyez patient.

À la fin de l'alignement, Bowtie2 renvoie des informations qui ressemblent à :

```
6757072 reads; of these:
  6757072 (100.00%) were unpaired; of these:
    1129248 (16.71%) aligned 0 times
    5164196 (76.43%) aligned exactly 1 time
    463628 (6.86%) aligned >1 times
83.29% overall alignment rate
```
On obtient ainsi :
- le nombre total de *reads* lus dans le fichier `.fastq.gz` (ici `6757072`)
- le nombre de *reads* non alignés « *aligned 0 times* » (`1129248`, soit `16.71%` du nombre total de *reads*)
- le nombre de *reads* alignés une seule fois (`5164196`)
- le nombre de *reads* alignés plus d'une fois (`463628`)
- un taux d'alignement global (`83.29%`)

Il faut être prudent si le nombre de *reads* non alignés est trop important (> 20%).


### Conversion des *reads* alignés en binaire, tri et indexation

Vous allez maintenant utiliser SAMtools pour :

1. Convertir le fichier `.sam` (fichier texte) créé par Bowtie2 en fichier `.bam`, qui est un fichier binaire, et qui donc prend moins de place sur le disque.
    ```
    $ samtools view -@ 2 -b map/bowtie.sam > map/bowtie.bam
    ```
    Cette étape va prendre plusieurs minutes. Comme votre machine dispose de 4 coeurs, nous allons en utiliser 2 (`-@ 2`) pour accélérer le calcul.  
    Comparez la taille deux fichiers `map/bowtie.sam` et `map/bowtie.bam`. Quel est le ratio de compression entre les deux formats de fichiers ?

2. Trier les *reads* alignés suivant l'ordre dans lequel ils apparaissent dans le génome.
    ```
    $ samtools sort -@ 2 map/bowtie.bam -o map/bowtie.sorted.bam
    ```
    Cette étape va prendre plusieurs minutes.

3. Indexer le fichier `.bam`. Cette étape est **indispensable** pour visualiser l'alignement avec IGV.
    ```
    $ samtools index map/bowtie.sorted.bam
    ```


### Comptage des *reads* alignés sur les gènes de *O. tauri*

Le comptage des *reads* alignés sur les gènes se fait avec HTSeq.

Toujours depuis votre shell Ubuntu et dans le répertoire `/mnt/c/Users/omics/rnaseq_tauri`, créez le répertoire `count` :
```
$ mkdir -p count
```

Puis lancez la commande (en une seule ligne) pour compter les *reads* alignés :
```
$ htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam map/bowtie.sorted.bam genome/GCF_000214015.3_version_140606.gff > count/count.txt
```

HTSeq renvoie le nombre d'annotations trouvées dans le fichier `.gff` puis affiche une progression de l'analyse. Les options du programme `htseq-count` sont décrites dans la [documentation](http://gensoft.pasteur.fr/docs/HTSeq/0.9.1/count.html).

Déterminez le nombre de *reads* alignés sur le gène `ostta18g01980`. Pour cela, lancez la commande :
```
$ grep ostta18g01980 count/count.txt
```
Vous pouvez aussi ouvrir le fichier `count/count.txt` avec la commande `less` puis chercher le gène `ostta18g01980` en tapant `/ostta18g01980` puis la touche <kbd>Entrée</kbd> (et enfin la touche <kbd>Q</kbd> pour quitter).


### Visualisation des *reads* alignés avec IGV

Pour visualiser l'alignement des *reads* sur le génome de référence avec IGV, vous avez besoin des fichiers suivants :
- Le génome de référence (`genome/GCF_000214015.3_version_140606.fna`).
- Les annotations du génome de référence (`genome/GCF_000214015.3_version_140606.gff`).
- Le fichier bam trié (`map/bowtie.sorted.bam`).
- L'index du bam trié (`map/bowtie.sorted.bam.bai`).

Lancez IGV et visualisez l'alignement des *reads* sur le génome de référence. Si vous avez oublié comme faire, visionnez la vidéo sur ce sujet qui vous a été proposée précédemment.

Visualisez particulièrement le gène `ostta18g01980`.


## 3.3 Automatisation de l'analyse : niveau 1

Tout cela est très bien mais les fichiers que vous avez générés (`map/bowtie.bam`, `map/bowtie.sorted.bam`, `cout/count.txt`...) portent des noms qui ne sont pas très informatifs sur l'échantillon dont ils proviennent.

Par ailleurs, entrer toutes ces commandes à la main, les unes après les autres, est pénible et source d'erreurs. Et il y a fort à parier que vous aurez complètement oublié ces commandes dans 1 semaine, voire dans 1 heure. 🤯 C'est parfaitement normal, il n'y a absolument aucun intérêt à se souvenir de toutes ces commandes.

Pour répondre à ces deux problèmes, de gestion de données et d'automatisation, nous allons introduire les notions Bash de variables et de scripts.

Mais d'abord, faites un peu de ménage en supprimant les fichiers créés précédemment :
```
$ rm -f reads_qc/*fastqc* index/*bt2 map/bowtie* count/count*
```

💣 Attention à l'utilisation de la commande `rm` qui supprime définitivement les fichiers.


### Variables

Une variable va simplement contenir de l'information qui sera utilisable autant de fois que nécessaire.

Création de variables :

```bash
$ toto=33
$ t="salut"
```

Attention : Il faut coller le nom de la variable et son contenu au symbole `=`.

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

### Script

Un script est un fichier texte qui contient des instructions Bash. Par convention, il porte l'extension `.sh`. L'objectif premier d'un script Bash est d'automatisr l'exécution de plusieurs commandes Bash, la plupart du temps pour manipuler ou analyser des fichiers.

Dans un script Bash, tout ce qui suit le symbole `#` est considéré comme un commentaire et n'est donc pas traité par Bash.


### Analyse RNA-seq

Testez le script `script1.sh` sur **un seul** de vos échantillons. Pour cela :

- Téléchargez le script `script1.sh` dans votre répertoire `rnaseq_tauri` avec la commande :
    ```
    $ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script1.sh
    ```

- Ouvrez le script `script1.sh` avec `nano`. Essayez de comprendre son fonctionnement, notamment l'utilisation des variables.

- Sur la deuxième ligne, modifiez la variable `sample` avec votre numéro d'échantillon. Sauvegardez le script (`ctrl + o`) et quittez nano (`ctrl + x`).  
    Rappel : pas d'espace avant et après le symbole `=` !

- Lancez le script avec la commande :
    ```
    $ bash script1.sh
    ```

Vérifiez que le déroulement du script se passe bien. Vous avez le temps de prendre un café (~ 20 ') ☕. Voir plusieurs ☕ 🍪 ☕ 🍪.

Évaluez approximativement le temps nécessaire au script 1 pour s'exécuter. ⏱️ À partir de cette valeur, extrapoler le temps nécessaire qu'il faudrait pour analyser les 3 échantillons.

Utilisez enfin la commande `tree` pour contempler votre travail (ici avec l'échantillon SRR2960338) :
```
$ tree
.
├── count
│   └── count-3.txt
├── genome
│   ├── GCF_000214015.3_version_140606.fna
│   └── GCF_000214015.3_version_140606.gff
├── index
│   ├── O_tauri.1.bt2
│   ├── O_tauri.2.bt2
│   ├── O_tauri.3.bt2
│   ├── O_tauri.4.bt2
│   ├── O_tauri.rev.1.bt2
│   └── O_tauri.rev.2.bt2
├── map
│   ├── bowtie-3.sorted.bam
│   └── bowtie-3.sorted.bam.bai
├── md5sum.txt
├── reads
│   ├── HCA-3_R1.fastq.gz
│   ├── HCA-3_R1_fastqc.html
│   ├── HCA-3_R1_fastqc.zip
│   ├── HCA-4_R1.fastq.gz
│   └── HCA-5_R1.fastq.gz
└── script1.sh

5 directories, 18 files
```


## 3.4 Automatisation de l'analyse : niveau 2

Le script précédent est pratique mais il ne conserve pas les informations liées à l'alignement générées par Bowtie2 (nombre de *reads* non-alignés, alignés une fois...).

Le [script 2](script2.sh) répond à ce problème. Téléchargez-le avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script2.sh
```

Ouvrez ce script avec `nano`. Vous remarquerez que la solution proposée pour conserver les informations liées à l'alignement est un peu particulière. Nous allons en discuter, mais dans un premier temps essayer de comprendre l'explication donnée [ici](https://stackoverflow.com/questions/876239/how-can-i-redirect-and-append-both-stdout-and-stderr-to-a-file-with-bash).


## 3.5 Automatisation de l'analyse : niveau 3 (ninja)

Le script précédent est intéressant mais il ne prend en compte qu'un seul échantillon à la fois. Quel ennui !

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

Le script 3 utilise une boucle pour automatiser l'analyse de plusieurs échantillons. Téléchargez-le avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script3.sh
```

Ouvrez ce script avec `nano`. Observez la structure du script et essayez de comprendre son fonctionnement.

La ligne `set -euo pipefail` tout au début du script va arrêter celui-ci :
- à la première erreur ;
- si une variable n'est pas définie ;
- si une erreur est rencontrée dans une commande avec un pipe (`|`).

C'est une mesure de sécurité importante pour votre script. Si vous le souhaitez, vous pouvez lire l'article de Aaron Maxwell à ce sujet : [Use the Unofficial Bash Strict Mode (Unless You Looove Debugging)](http://redsymbol.net/articles/unofficial-bash-strict-mode/)


Toujours avec `nano`, modifiez le script 3 avec les numéros d'échantillons que vous avez à analyser. Faites bien attention à la variable concernée et à sa syntaxe.

Si vous pensez en avoir le temps, lancez le script 3. Comme ce script va automatiser toute l'analyse, il va fonctionner plusieurs dizaines de minutes et vous aurez peut-être besoin de fermez votre terminal. Pour ne pas arrêter brutalement l'analyse, lancez le script de cette manière :

```
$ nohup bash script3.sh &
```
Mais pour autant, n'éteignez pas votre ordinateur !


Le message 
```
nohup: ignoring input and appending output to 'nohup.out'
```
vous rappelle que les messages qui apparaissaient habituellement à l'écran seront redirigés dans le fichier `nohup.out`.


##  3.6 Comparaison avec les logiciels utilisés dans Galaxy (si vous avez du temps)

Connectez-vous maintenant à votre compte sur Galaxy. Essayez de retrouver les versions des logiciels que vous utilisés (FastQC, Bowtie2, SAMtools, HTSeq).

Pour ce faire, dans votre *History*, cliquez sur le nom d'un résultat d'analyse, puis cliquez sur le petit i entouré (ℹ️) et lisez les informations de la section *Job Dependencies*.

Comparez les versions des logiciels disponibles dans Galaxy et de ceux que vous avez utilisé sur votre machine.
