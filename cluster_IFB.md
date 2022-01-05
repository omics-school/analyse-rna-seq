---
title: Calcul sur le cluster de l'IFB
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---


Dans cette activité, vous allez analyser les données RNA-seq de *O. tauri* avec le cluster *National Network of Computational Resources* (NNCR) de l'Institut Français de Bioinformatique (IFB). Ce cluster utilise un système d'exploitation Linux.


## Remarques préables

L'accès au cluster de l'IFB vous est fourni dans le cadre du DU Omiques. Cet accès sera révoqué à l'issue de la formation. 

Si à l'issue de cette formation, vous souhaitez continuer à utiliser ce cluster pour votre projet de recherche, connectez-vous sur votre [interface](https://my.cluster.france-bioinformatique.fr/manager2/project) puis cliquez sur le bouton *Request A New Project* et précisez en quelques mots votre projet. Plusieurs utilisateurs peuvent être associées à un même projet et partager des données.

Si vous avez besoin d'un logiciel spécifique sur le cluster. N'hésitez pas à le demander gentillement sur le site [Cluster Community Support](https://community.cluster.france-bioinformatique.fr/). Les administrateurs sont en général très réactifs.


## 0. Connexion au cluster

Sous Windows, ouvrez un terminal Ubuntu. Si vous avez oublié comment faire, consultez le [tutoriel Unix](https://omics-school.github.io/unix-tutorial/tutoriel/README).

Déplacez vous ensuite dans le répertoire `/mnt/c/Users/omics` :

```bash
$ cd /mnt/c/Users/omics
```

🔔 Rappels : 

- Ne tapez pas le `$` en début de ligne et faites attention aux majuscules et aux minuscules (surtout pour `Users`) !
- Utilisez le copier / coller.
- Utilisez la complétion des noms de fichier et de répertoires avec la touche <kbd>Tab</kbd>.

Connectez-vous en SSH au cluster avec les identifiants (login et mot de passe) que vous avez du recevoir par e-mail.

La syntaxe est de la forme :
```bash
$ ssh LOGIN@core.cluster.france-bioinformatique.fr
```

avec `LOGIN` votre identifiant sur le cluster. 

Si c'est la première fois que vous vous connectez au cluster, répondez `yes` à la question 
```
Are you sure you want to continue connecting (yes/no)?
```

Entrez ensuite votre mot de passe en **aveugle**, c'est-à-dire sans qu'aucun caractère ne soit affiché à l'écran. C'est assez déstabilisant la première fois puis on s'habitue.

🔔 **Attention** 🔔 Le cluster est protégé contre certaines attaques. Si vous entrez un mot de passe erronné plusieurs fois de suite, votre IP va être bannie et vous ne pourrez plus vous connecter (temporairement) au serveur.


Pour vous déconnecter du cluster et revenir à votre terminal local, pressez la combinaison de touches <kbd>Ctrl</kbd>+<kbd>D</kbd>.

Un cluster est un ensemble de machines. La machine à laquelle vous venez de vous connecter est le noeud de connexion. C'est aussi depuis cette machine que vous lancerez vos analyses. 

**Remarque :** Vous lancerez vos calculs **depuis** le noeud de connexion mais pas **sur** le noeud de connexion. Il est interdit de lancer une analyse sur le noeud de connexion sous peine de voir votre compte suspendu.


## 1. Environnement logiciel 

Si vous vous êtes déconnectés du cluster, reconnectez-vous avec la commande `ssh` précédente.

Par défaut, aucun logiciel de bioinformatique n'est présent. Pour vous en convaincre, essayez de lancer la commande :
```bash
$ bowtie2 --version
```
Vous devriez obtenir un message d'erreur du type : `-bash: bowtie2 : commande introuvable`

Chaque logiciel doit donc être chargé individuellement avec l'outil `module`.

Utilisez la commande suivante pour compter le nombre de logiciels disponibles avec `module` :
```bash
$ module avail -l | wc -l
```

Chargez ensuite les logiciels `fastqc`, `bowtie2`, `samtools` et `htseq` avec les commandes suivantes :
```bash
$ module load fastqc/0.11.9
$ module load bowtie2/2.3.5
$ module load samtools/1.9
$ module load htseq/0.11.3
```

Vérifiez que les logiciels sont bien disponibles en affichant leurs versions :

```bash
$ fastqc --version
FastQC v0.11.9
```

```bash
$ bowtie2 --version
/home/duo/miniconda3/envs/rnaseq/bin/bowtie2-align-s version 2.3.5.1
64-bit
Built on
Wed Apr 17 02:40:25 UTC 2019
[...]
```

```bash
$ samtools --version
samtools 1.9
Using htslib 1.9
Copyright (C) 2018 Genome Research Ltd.
```

```bash
$ htseq-count -h | tail -n 4
Written by Simon Anders (sanders@fs.tum.de), European Molecular Biology
Laboratory (EMBL) and Fabio Zanini (fabio.zanini@stanford.edu), Stanford
University. (c) 2010-2019. Released under the terms of the GNU General Public
License v3. Part of the 'HTSeq' framework, version 0.11.3.
```

## 2. Stockage des données

Votre répertoire utilisateur sur le noeud de connexion (`/shared/home/LOGIN`) ne doit pas contenir de donnée car l'espace disponible est limité à 100 Go. Un espace de stockage a été créé pour vous dans le répertoire  `/shared/projects/form_2021_29/LOGIN` (avec `LOGIN` votre identifiant sur le cluster). Par la suite, cet espace sera appelé « répertoire de travail ».

De plus, le répertoire `/shared/projects/form_2021_29/data/rnaseq_tauri` contient les données dont vous aurez besoin pour ce projet. Vous n'avez accès à ce répertoire qu'en lecture, c'est-à-dire que vous pouvez seulement parcourir les répertoires et lire les fichiers de ce répertoire (pas de modification, d'ajout ou de suppression).

De quels fichiers avez-vous besoin pour l'analyse des données RNA-seq de *O. tauri* ? 

Vérifiez que tous les fichiers nécessaires sont bien présents dans `/shared/projects/form_2021_29/data/rnaseq_tauri`.

Vérifiez l'intégrité des fichiers `.fastq.gz` situés dans le répertoire `/shared/projects/form_2021_29/data/rnaseq_tauri/reads` avec les commandes suivantes :

```bash
$ cd /shared/projects/form_2021_29/data/rnaseq_tauri
```

*Rappel : n'entrez pas le symbole $ en début de ligne*

puis 

```bash
$ srun -A form_2021_29 md5sum -c reads_md5sum.txt
```

N'oubliez pas le `srun -A form_2021_29` en début de commande :

- L'instruction `srun` est spécifique au cluster. 
- L'option `-A form_2021_29` spécifie quel projet utiliser (facturer) pour cette commande. Un même utilisateur peut appartenir à plusieurs projets. Le nombre d'heures de calcul attribuées à un projet étant limité, il est important de savoir quel projet imputé pour telle ou telle commande. Pensez-y pour vos futurs projets.


Déplacez-vous maintenant dans votre répertoire de travail `/shared/projects/form_2021_29/LOGIN` (avec `LOGIN` votre identifiant sur le cluster).

Créez le répertoire `rnaseq_tauri` et déplacez-vous à l'intérieur. Dorénavant vous ne travaillerez plus qu'à partir de ce répertoire.

La commande `pwd` devrait vous renvoyer quelque chose du type :

```bash
$ pwd
/shared/projects/form_2021_29/LOGIN/rnaseq_tauri`
```

avec `LOGIN` votre identifiant sur le cluster. 🆘 Appelez à l'aide si vous ne parvenez pas à être dans le bon répertoire.


## 3.1 Analyse d'un échantillon

**Remarque préalable** : l'indexation du génome de référence avec le logiciel `bowtie2` a déjà été effectué pour vous. Pour vous en convraincre, affichez le contenu du répertoire `/shared/projects/form_2021_29/data/rnaseq_tauri/genome` et vérifiez l'existence de fichiers avec l'extension `.bt2`, spécifiques des fichiers index créés par `bowtie2`.

Depuis le cluster de l'IFB, vérifiez que vous êtes toujours dans votre répertoire `/shared/projects/form_2021_29/LOGIN/rnaseq_tauri`.

Téléchargez ensuite le script `script4.sh` avec la commande :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script4.sh
```

Lancez ce script avec la commande :

```bash
$ sbatch -A form_2021_29 script4.sh
```

Vous devriez obtenir un message du type `Submitted batch job 20716345`. Ici, `20716345` correspond au numéro du job.

Notez bien le numéro de votre job.

Vérifiez que votre script est en train de tourner avec la commande :

```bash
$ squeue -u $USER
```

**Remarque** Voici quelques statuts (colonne `ST`) de job intéressant :

- `CA` (*cancelled*) : le job a été annulé
- `F` (*failled*) : le job a planté
- `PD` (*pending*) : le job est en attente que des ressources soient disponibles
- `R` (*running*) : le job est lancé


Et pour avoir plus de détails, utilisez la commande :
```bash
$ sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j JOBID
```
avec `JOBID` (en fin de ligne) le numéro de votre job à remplacer par le vôtre.

Voici un exemple de sortie que vous pourriez obtenir :

```bash
$ sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j 20716345
       JobID    JobName      State               Start    Elapsed    CPUTime        NodeList 
------------ ---------- ---------- ------------------- ---------- ---------- --------------- 
20716345     script4.sh    RUNNING 2022-01-05T23:37:36   00:08:58   00:08:58     cpu-node-24 
20716345.ba+      batch    RUNNING 2022-01-05T23:37:36   00:08:58   00:08:58     cpu-node-24 
20716345.0       fastqc    RUNNING 2022-01-05T23:37:37   00:00:54   00:00:34     cpu-node-24 
```

Si vous affichez le contenu de votre répertoire courant, vous devriez voir l'apparition d'un fichier `slurm-JOBID.out` où `JOBID` est le numéro de votre job. Ce fichier contient la sortie, c'est-à-dire le *log* de votre script.

Affichez son contenu avec la commande `cat`. Par exemple :

```bash
$ cat slurm-JOBID.out
```

avec `JOBID` le numéro de votre job.

Si vous attendez 1 ou 2 minutes et relancez la commande `sacct` précédente, votre job a du passer à une nouvelle étape.

Par exemple :

```bash
$ sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j 20716345
       JobID    JobName      State               Start    Elapsed    CPUTime        NodeList 
------------ ---------- ---------- ------------------- ---------- ---------- --------------- 
20716345     script4.sh    RUNNING 2022-01-05T23:37:36   00:05:22   00:05:22     cpu-node-24 
20716345.ba+      batch    RUNNING 2022-01-05T23:37:36   00:05:22   00:05:22     cpu-node-24 
20716345.0       fastqc  COMPLETED 2022-01-05T23:37:37   00:00:54   00:00:54     cpu-node-24 
20716345.1      bowtie2    RUNNING 2022-01-05T23:38:31   00:04:27   00:04:27     cpu-node-24 
```

Nous allons maintenant améliorer le script d'analyse. Annulez votre job en cours avec la commande :
```bash
$ scancel JOBID
```

où `JOBID` est le numéro de votre job.

Faites aussi un peu de ménage en supprimant les fichiers créés précédemment avec la commande :
```bash
$ rm -rf map/ reads_qc/ count/ slurm*.out
```

## 3.2 Analyse plus rapide d'un échantillon

L'objectif est maintenant « d'aller plus vite » en attribuant plusieurs coeurs pour l'étape d'alignement des reads sur le génome avec `bowtie2`.

Toujours depuis le cluster de l'IFB, dans le répertoire `rnaseq_tauri` de votre répertoire de travail, téléchargez le script `script5.sh` avec la commande :
```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script5.sh
```

Identifiez les différences avec le script précédent, par exemple avec la commande `diff` : 
```bash
$ diff script4.sh script5.sh
```

Lancez ensuite votre analyse :
```bash
$ sbatch -A form_2021_29 script5.sh
```

Notez bien le numéro de job renvoyé.

Vérifiez que votre job est bien lancé avec la commande :
```bash
$ squeue -u $USER
```

Le fichier `slurm-JOBID.out` est également créé et contient les sorties du script. Pour consulter son contenu, tapez :
```bash
$ cat slurm-JOBID.out
```

avec `JOBID` le numéro de votre job.


Suivez également en temps réel l'exécution de votre job avec la commande :
```bash
$ watch sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j JOBID
```
avec `JOBID` le numéro de votre job.

Remarques : 

- La commande `watch` est utilisée ici pour « surveiller » le résultat de la commande `sacct`.
- L'affichage est rafraichi toutes les 2 secondes.

Votre job devrait prendre une petite dizaine de minutes pour se terminer. Laissez le cluster travailler et profitez-en pour vous préparer un thé ou un café.

Quand les status (colonne `State`) du job et de tous les job steps sont à `COMPLETED`, stoppez la commande `watch` en appuyant sur la combinaison de touches <kbd>Ctrl</kbd> + <kbd>C</kbd>.

Vérifiez avec la commande `tree` que les fichiers suivants ont bien été créés :


- `HCA-37_R1_fastqc.html`
- `HCA-37_R1_fastqc.zip`
- `bowtie-37.sorted.bam`
- `count-37.txt`
- `slurm-jobID.out` (avec `jobID` le numéro de votre job)

Vérifiez que la somme de contrôle du fichier `count-37.txt` est bien `cbc9ff7ed002813e16093332c7abfed4`.

## 3.3 Analyse de plusieurs échantillons

Toujours depuis le cluster de l'IFB, dans le répertoire `rnaseq` de votre répertoire de travail, téléchargez le script `script6.sh` avec la commande :
```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script6.sh
```

Nous pourrions analyser d'un seul coup les 47 échantillons (fichiers `.fastq.gz`) mais pour ne pas consommer trop de ressources sur le cluster, nous allons limiter notre analyse à 4 échantillons.

Lancez votre analyse avec la commande :
```bash
$ sbatch script6.sh
```

Notez bien le numéro de job renvoyé.

Vous pouvez suivre en temps réel l'exécution de votre job avec la commande :
```bash
$ watch sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j jobID
```
avec `jobID` le numéro de votre job.

Patientez une dizaine de minutes que tous les jobs et job steps soient terminés. 

Quand les status (colonne `State`) de tous les jobs et job steps sont à `COMPLETED`, stoppez la commande `watch` en appuyant sur la combinaison de touches <kbd>Ctrl</kbd> + <kbd>C</kbd>.


## 4. L'heure de faire les comptes

Expérimentez la commande `sreport` pour avoir une idée du temps de calcul consommé par tous vos jobs :

```bash
$ sreport -t hour Cluster UserUtilizationByAccount Start=2020-01-01 End=$(date --iso-8601)T23:59:59 Users=$USER
```

La colonne `Used` indique le nombre d'heures de temps CPU consommées. Cette valeur est utile pour estimer le « coût CPU » d'un projet.

Voici un exemple de rapport produit par `sreport` :

```bash
$ sreport -t hour Cluster UserUtilizationByAccount Start=2020-01-01 End=$(date --iso-8601)T23:59:59 Users=$USER
--------------------------------------------------------------------------------
Cluster/User/Account Utilization 2020-01-01T00:00:00 - 2021-03-18T21:59:59 (38268000 secs)
Usage reported in CPU Hours
--------------------------------------------------------------------------------
  Cluster     Login     Proper Name         Account     Used   Energy 
--------- --------- --------------- --------------- -------- -------- 
     core  ppoulain  Pierre Poulain       dubii2021      400        0 
     core  ppoulain  Pierre Poulain     du_bii_2019      246        0 
     core  ppoulain  Pierre Poulain uparis_duo_2020      129        0 
     core  ppoulain  Pierre Poulain        minomics        5        0 
```

Ainsi, l'utilisateur `ppoulain` a déjà consommé 129 heures de temps CPU sur le projet `uparis_duo_2020`.


## 5. Récupération des données

### 5.1 scp

⚠️ Pour récupérer des fichiers sur le cluster en ligne de commande, vous devez lancer la commande `scp` depuis un shell Unix sur votre machine locale. ⚠️

Depuis un shell Unix sur votre machine locale, déplacez-vous dans le répertoire `/mnt/c/Users/omics` et créez le répertoire `rnaseq_cluster`. 

Déplacez-vous dans ce nouveau répertoire.

Utilisez la commande `pwd` pour vérifier que vous êtes bien dans le répertoire `/mnt/c/Users/omics/rnaseq_cluster`. 

Lancez ensuite la commande suivante pour récupérer les fichiers de comptage :
```bash
$ scp login@core.cluster.france-bioinformatique.fr:/shared/projects/uparis_duo_2020/login/rnaseq/count*.txt .
```

où `login` est votre identifiant sur le cluster. Faites bien attention à garder le `.` tout à la fin de la commande.


Vérifiez que la somme de contrôle MD5 du fichier `count-37.txt` est bien la même que précédemment.


### 5.2 Filezilla

Lancez le logiciel FileZilla ([comme ceci](img/filezilla.png)). Puis entrez les informations suivantes :

- Hôte : `sftp://core.cluster.france-bioinformatique.fr`
- Identifiant : votre login sur le cluster
- Mot de passe : votre mot de passe sur le cluster

Cliquez ensuite sur le bouton *Connexion rapide*. Cliquez sur *OK* dans la fenêtre *Clé de l'hôte inconnue*

Une fois connecté, dans le champs texte à coté de *Site distant* (à droite de la fenêtre), entrez le chemin `/shared/projects/uparis_duo_2020/` voire directement votre répertoire de travail `/shared/projects/uparis_duo_2020/login` (avec `login` votre identifiant sur le cluster).

Essayez de transférer des fichiers dans un sens puis dans l'autre. Double-cliquez sur les fichiers pour lancer les transferts.



