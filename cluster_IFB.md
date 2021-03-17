---
title: Calcul sur le cluster de l'IFB
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---


Dans cette activité, vous allez analyser les données RNA-seq de *O. tauri* avec le cluster *National Network of Computational Resources* (NNCR) de l'Institut Français de Bioinformatique (IFB). Ce cluster utilise un système d'exploitation Linux.


## Remarques préables

L'accès au cluster de l'IFB vous est fourni dans le cadre du DU Omiques. Cet accès sera révoqué à l'issue de la formation. 

Si vous souhaitez continuer à utiliser ce cluster pour votre projet, connectez-vous sur votre [interface](https://my.cluster.france-bioinformatique.fr/manager2/project) puis cliquez sur le bouton *Request A New Project* et précisez en quelques mots votre projet. Plusieurs utilisateurs peuvent être associées à un même projet et partager des données.

Si vous avez besoin d'un logiciel spécifique sur le cluster. N'hésitez pas à le demander sur le site [Cluster Community Support](https://community.cluster.france-bioinformatique.fr/). Les administrateurs sont en général très réactifs.


## 0. Connexion au cluster

Depuis votre station de travail, ouvrez un shell Unix.

Connectez-vous en SSH au cluster avec les identifiants (login et mot de passe) que vous avez du recevoir par e-mail.

La syntaxe est de la forme :
```
ssh login@core.cluster.france-bioinformatique.fr
```

avec `login` votre identifiant. 

Si c'est la première fois que vous vous connectez au cluster, répondez `yes` à la question 
```
Are you sure you want to continue connecting (yes/no)?
```

Vous entrerez ensuite votre mot de passe en aveugle, c'est-à-dire qu'aucun caractère ne sera affiché à l'écran. C'est assez déstabilisant la première fois puis on s'habitue.

Un cluster est un ensemble de machines. La machine à laquelle vous venez de vous connecter est le noeud de connexion. C'est aussi depuis cette machine que vous lancerez vos analyses. 

**Remarque :** Vous lancerez vos calculs **depuis** le noeud de connexion mais pas **sur** le noeud de connexion. Il est interdit de lancer une analyse sur le noeud de connexion sous peine de voir votre compte suspendu.


## 1. Stockage des données

Votre répertoire utilisateur sur le noeud de connexion (`/shared/home/login`) ne doit pas contenir vos données car l'espace disponible est limité à 100 Go. Un espace de stockage a été créé pour vous dans le répertoire  `/shared/projects/uparis_duo_2020/login` (avec `login` votre identifiant sur le cluster). Par la suite, cet espace sera appelé « répertoire de travail ».

De plus, le répertoire `/shared/projects/uparis_duo_2020/data` contient les données dont vous aurez besoin pour ce projet. Vous n'avez accès à ce répertoire qu'en lecture, c'est-à-dire que vous pouvez seulement parcourir les répertoires et lire les fichiers de ce répertoire (pas de modification, d'ajout ou de suppression).

De quels fichiers aviez-vous besoin pour l'analyse des données RNA-seq de *O. tauri* ? 

Vérifiez que tous les fichiers nécessaires sont bien présents dans `/shared/projects/uparis_duo_2020/data`.

Vérifiez l'intégrité des fichiers `.fastq.gz` situés dans le répertoire `/shared/projects/uparis_duo_2020/data/reads` avec les commandes suivantes :

```
$ cd /shared/projects/uparis_duo_2020/data/reads
```

*Rappel : n'entrez pas le symbole $ en début de ligne*

puis 

```
$ srun md5sum -c md5sum.txt
```

N'oubliez pas le `srun` en début de commande, sans quoi vous allez recevoir un appel énervé de l'administrateur du cluster.


Déplacez-vous maintenant dans votre répertoire de travail `/shared/projects/uparis_duo_2020/login` (avec `login` votre identifiant sur le cluster).

Créez le répertoire `rnaseq` et déplacez-vous à l'intérieur. Dorénavant vous ne travaillerez qu'à partir de ce répertoire.


## 2. Environnement logiciel 

Par défaut, aucun logiciel de bioinformatique n'est présent. Pour vous en convaincre, essayez de lancer la commande :
```
$ bowtie2 --version
```
Vous devriez obtenir un message d'erreur du type :
```
-bash: bowtie2 : commande introuvable
```

Chaque logiciel doit donc être chargé individuellement avec l'outil `module`.

Utilisez la commande suivante pour compter le nombre de logiciels disponibles avec `module` :
```
$ module avail -1 | wc -l
```

Chargez ensuite les logiciels `fastqc`, `bowtie2`, `samtools` et `htseq` avec les commandes suivantes :
```
$ module load fastqc/0.11.9
$ module load bowtie2/2.3.5
$ module load samtools/1.9
$ module load htseq/0.11.3
```

Vérifiez que les logiciels sont bien installés en affichant leurs versions :

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
$ htseq-count -h | tail -4
Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology
Laboratory (EMBL) and Fabio Zanini (fabio.zanini@stanford.edu), Stanford
University. (c) 2010-2019. Released under the terms of the GNU General Public
License v3. Part of the 'HTSeq' framework, version 0.11.3.
```

## 3.1 Analyse d'un échantillon

Depuis le cluster de l'IFB, dans le répertoire `rnaseq` de votre répertoire de travail, téléchargez le script `script4.sh` avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script4.sh
```

Lancez ensuite ce script avec la commande :
```
$ sbatch script4.sh
```

Notez bien le numéro de job renvoyé.

Vérifiez que votre script est en train de tourner avec la commande :

```
$ squeue -u $USER
```

**Remarque** Voici quelques statuts (colonne `ST`) de job intéressant :

- `CA` (*cancelled*) : le job a été annulé
- `F` (*failled*) : le job a planté
- `PD` (*pending*) : le job est en attente que des ressources soient disponibles
- `R` (*running*) : le job est lancé


Et pour avoir plus de détails :
```
$ sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j jobID
```
avec `jobID` le numéro de votre job.


Annulez votre job avec la commande :
```
$ scancel jobID
```

où `jobID` est le numéro de votre job.

Faites un peu de ménage en supprimant les fichiers créés avec la commande :
```
$ rm -f bowtie*bam HCA*html HCA*zip count*txt
```

## 3.2 Analyse d'un échantillon plus rapide

L'objectif est maintenant « d'aller plus vite » en attribuant plusieurs coeurs pour l'étape d'alignement des reads sur le génome avec `bowtie2`.

Toujours depuis le cluster de l'IFB, dans le répertoire `rnaseq` de votre répertoire de travail, téléchargez le script `script5.sh` avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script5.sh
```

Identifiez les différences avec le script précédent.

Lancez ensuite votre analyse :
```
$ sbatch script5.sh
```

Notez bien le numéro de job renvoyé.

Vérifiez que votre job est bien lancé avec la commande :
```
$ squeue -u $USER
```

Le fichier `slurm-jobID.out` est également créé et contient les sorties du script. Pour le consultez en temps réél, tapez :
```
$ tail -f slurm-jobID.out
```

avec `jobID` le numéro de votre job.

Pour quitter, appuyez sur la combinaison de touches <kbd>Ctrl</kbd> + <kbd>C</kbd>.


Suivez en temps réel l'exécution de votre job avec la commande :
```
$ watch sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j jobID
```
avec `jobID` le numéro de votre job.

Remarques : 

- La commande `watch` est utilisée ici pour « surveiller » le résultat de la commande `sacct`.
- L'affichage est rafraichi toutes les 2 secondes.

Votre job devrait prendre une petite dizaine de minutes pour se terminer. Laissez le cluster travailler et profitez-en pour vous faire un thé ou un café.

Quand les status (colonne `State`) du job et de tous les job steps sont à `COMPLETED`, stoppez la commande `watch` en appuyant sur la combinaison de touches <kbd>Ctrl</kbd> + <kbd>C</kbd>.

Vérifiez que les fichiers suivants ont bien été créés dans votre répertoire :

- `HCA-37_R1_fastqc.html`
- `HCA-37_R1_fastqc.zip`
- `bowtie-37.sorted.bam`
- `count-37.txt`
- `slurm-jobID.out` (avec `jobID` le numéro de votre job)

Vérifiez que la somme de contrôle du fichier `count-37.txt` est bien `cbc9ff7ed002813e16093332c7abfed4`.

## 3.3 Analyse de plusieurs échantillons

Toujours depuis le cluster de l'IFB, dans le répertoire `rnaseq` de votre répertoire de travail, téléchargez le script `script6.sh` avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script6.sh
```

Nous pourrions analyser d'un seul coup les 47 échantillons (fichiers `.fastq.gz`) mais pour ne pas consommer de trop de ressources sur le cluster, nous allons limiter notre analyse à 4 échantillons.

Lancez votre analyse avec la commande :
```
$ sbatch script6.sh
```

Notez bien le numéro de job renvoyé.

Vous pouvez suivre en temps réel l'exécution de votre job avec la commande :
```
$ watch sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j jobID
```
avec `jobID` le numéro de votre job.

Patientez une dizaine de minutes que tous les jobs et job steps soient terminées. 

Quand les status (colonne `State`) du job et de tous les job steps sont à `COMPLETED`, stoppez la commande `watch` en appuyant sur la combinaison de touches <kbd>Ctrl</kbd> + <kbd>C</kbd>.

## 4. L'heure de faire les comptes

Expérimentez la commande `sreport` pour avoir une idée du temps de calcul consommé par tous vos jobs :

```
$ sreport Cluster UserUtilizationByAccount Start=2020-01-01 Users=$USER
```


## 5. Récupération des données

### 5.1 scp

⚠️ Pour récupérer des fichiers sur le cluster en ligne de commande, vous devez lancer la commande `scp` depuis un shell Unix sur votre machine locale. ⚠️

Depuis un shell Unix sur votre machine locale, déplacez-vous dans le répertoire `/mnt/c/Users/omics` et créez le répertoire `rnaseq_cluster`. 

Déplacez-vous dans ce nouveau répertoire.

Utilisez la commande `pwd` pour vérifier que vous êtes bien dans le répertoire `/mnt/c/Users/omics/rnaseq_cluster`. 

Lancez ensuite la commande suivante pour récupérer les fichiers de comptage :
```
$ scp login@core.cluster.france-bioinformatique.fr:/shared/projects/uparis_duo_2020/login/rnaseq/count*.txt .
```

où `login` est votre identifiant sur le cluster. Faites bien attention à garder le `.` tout à la fin de la commande.


Vérifiez que la somme de contrôle MD5 du fichier `count-37.txt` est bien le même que précédemment.


### 5.2 Filezilla

Lancez le logiciel FileZilla ([comme ceci](img/filezilla.png)).

Puis entrez les informations suivantes :

- Hôte : `sftp://clore.cluster.france-bioinformatique.fr`
- Identifiant : votre login sur le cluster
- Mot de passe : votre mot de passe sur le cluster

Cliquez ensuite sur le bouton *Connexion rapide*. Cliquez sur *OK* dans la fenêtre *Clé de l'hôte inconnue*

Une fois connecté, dans le champs texte à coté de *Site distant* (à droite de la fenêtre), entrez le chemin `/shared/projects/uparis_duo_2021/` voire directement votre répertoire de travail `/shared/projects/uparis_duo_2021/login` (avec `login` votre identifiant sur le cluster).

Essayez de transférer des fichiers dans un sens puis dans l'autre. Double-cliquez sur les fichiers pour lancer les transferts.



