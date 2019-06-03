---
title: Calcul sur le cluster de l'IFB
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---


Dans cette activité, vous allez analyser les données RNA-seq de *O. tauri* sur le cluster NNCR de l'Institut Français de Bioinformatique (IFB). Ce cluster utilise un système d'exploitation Unix / Linux.

# Remarques préables

L'accès au cluster de l'IFB vous est fourni dans le cadre du DU Omique. C'est accès sera révoqué à l'issue de la formation, fin janvier 2020. 

Si vous souhaitez continuer à utiliser ce cluster, faites-en la demande en remplissant le formulaire [IFB core cluster - account request](https://www.france-bioinformatique.fr/fr/ifb-core-cluster-account-request) en précisant en quelques mots votre projet.

Si vous avez besoin d'un logiciel spécifique sur le cluster. N'hésitez pas à le demander sur le site [Cluster Community Support](https://community.cluster.france-bioinformatique.fr/)


# Connexion 

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

Vous entrerez ensuite votre mot de passe en aveugle, c'est-à-dire qu'aucun caractère ne sera affiché à l'écran.


# Découverte de l'environnement

## Noeud de connexion

Un cluster est un ensemble de machines. La machine à laquelle vous venez de vous connecter est le noeud de connexion. C'est aussi depuis cette machine que vous lancerez vos analyses. 

**Remarque :** Vous lancerez vos calculs **depuis** le noeud de connexion mais pas **sur** la noeud de connexion. Il est interdit de lancer une analyse sur le noeud de connexion sous peine de voir votre compte suspendu.


## Stockage des données

Votre répertoire utilisateur sur le noeud de connexion (`/shared/home/login`) ne doit pas contenir vos données car l'espace disponible est limité à 100 Go. Un espace de stockage a été créé pour vous dans le répertoire  `/shared/projects/du_o_2019/login`. Par la suite, cet espace sera appelé « répertoire de travail ».

De plus, le répertoire `/shared/projects/du_o_2019/data` contient les données dont vous aurez besoin pour ce projet. Vous n'avez accès à ce répertoire qu'en lecture, c'est-à-dire que vous pouvez seulement parcourir les répertoires et lire les fichiers (pas de modification, d'ajout ou de suppression).

Vérifiez que tous les fichiers nécessaires pour l'analyse des données RNA-seq de *O. tauri* sont bien présents dans `/shared/projects/du_o_2019/data`.

Essayez de créer un fichier dans les répertoires :

-  `/shared/projects/du_o_2019/data`
-  `/shared/projects/du_o_2019/login`


## Environnement logiciel 

L'environnement logiciel nécessaire pour l'analyse RNA-seq a été installée par les administrateurs du cluster.

Pour l'activer, lancer la commande :
```
$ module ...
```

Vérifiez que `fastqc`, `bowtie2`, `samtools` et `htseq-count` sont disponibles. 

Quelles sont les versions de ces outils ? Si besoin, retournez voir le [Tutoriel de l'analyse RNA-seq](analyse_RNA-seq_O_tauri.md). 

Est-ce que sont les mêmes versions que sur le serveur du DU ?


## Préparation des données

Dans votre répertoire de travail (`/shared/projects/du_o_2019/login`), créez le répertoire `RNAseq`.

Copiez à l'intérieur de ce répertoire les fichiers dont vous aurez besoin pour travailler :

- le génome de référence
- les annotations du génome
- les 2 ou 3 fichiers de reads


## Commandes manuelles

Depuis le répertoire `RNAseq` de votre répertoire de travail, lancez un contrôle qualité d'un fichier de séquençage avec la commande :
```
$ srun fastqc nom-fichier-fastq.gz
```
où `nom-fichier-fastq.gz` est le fichier contenant l'échantillon que vous avez choisi.

La commande `srun` va lancer l'analyse du contrôle qualité (`fastqc nom-fichier-fastq.gz`) sur un des noeuds de calcul du cluster. 

FastQC va produire deux fichiers (un fichier avec l'extension `.html` et un autre avec l'extension `.zip`). Copiez le fichier `.html` sur votre machine locale avec le logiciel FileZilla ou la commande `scp`. Visualisez ce fichier avec votre navigateur web.


## Automatisation 1 

Depuis le cluster de l'IFB, dans le répertoire `RNAseq` de votre répertoire de travail, téléchargez le script 3 avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script3.sh
```

Remarquez que c'est exactement le même script qui fonctionnait sur le serveur du DU.

Pour que ce premier test soit assez rapide, modifiez-le avec `nano` pour l'adapter à **un seul** de vos échantillons.

Puis lancez-le avec la commande :
```
$ srun bash script3.sh
```

Pour lancer votre analyse puis fermer votre session (et partir en week-end), utilisez plutôt 
```
$ nohup srun bash script3.sh &
```

Pour vérifier l'état de votre job, appelez la commande :
```
$ squeue -u login
```

Par exemple :
```
             JOBID PARTITION     NAME     USER ST       TIME  NODES NODELIST(REASON)
            438536      fast     bash ppoulain  R       4:04      1 cpu-node-8
```

La colonne `ST` indique le statut de votre job. Si il est actif, son statut doit être `R` (pour *running*). La colonne `NODELIST(REASON)` indique sur quel noeud du cluster a été lancé votre job (ici `cpu-node-8`).

Par défaut, la commande `srun` va lancer votre job sur un noeud avec un seul CPU.

Si vous avez besoin de supprimer un de vos jobs, utilisez la commande 
```
$ scancel job-id
```

où `job-id` est l'identifiant de votre job (colonne `JOBID` indiquée par la commande `squeue`).

Supprimez un job que vous avez lancé.


@JULIEN : 

- avec la commande `nohup run ... &` un utilisateur peut lancer plusieurs jobs différents ?
- `srun` lancer le job immédiatement... si des ressources sont disponibles. Que se passe-t-il si ce n'est pas le cas ? Le job échoue ?



## Automatisation 2 (sbatch)

Toujours depuis le cluster de l'IFB, téléchargez le script 4 avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script4.sh
```

Quelles différences observez-vous avec le script précédent ?

Ouvrez ce fichier avec `nano` puis modifiez-le pour adapter votre adresse e-mail et vos numéros d'échantillons.

Lancez ensuite votre analyse :
```
$ sbatch script4.sh
```

Un message équivalent à `Submitted batch job 440893` vous indique que votre job a correctement été lancé et vous indiquer son numéro d'identification `440893`.

Vérifiez que votre job est bien lancé avec 
```
$ squeue -u login
```

Le fichier `slurm-440893.out` est également créé et contient les sorties du script. Pour le consultez en temps réél, tapez :
```
$ tail -f slurm-440893.out
```

Pour quitter, appuyez sur la combinaison de touches <kbd>Ctrl</kbd> + <kbd>C</kbd>.

@JULIEN :

- Dans un script sbatch, toutes les lignes *exécutables* doivent être lancées avec srun. Pourquoi ?


## Automatisation 3 (sbatch + multi-coeurs)

Toujours depuis le cluster de l'IFB, téléchargez le script 5 avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script5.sh
```

Quelles différences observez-vous avec le script précédent ?

Ouvrez ce fichier avec `nano` puis modifiez-le pour adapter votre adresse e-mail et vos numéros d'échantillons.

Lancez ensuite votre analyse :
```
$ sbatch script5.sh
```

Affichez en temps réel le fichier qui contient la sortie du script. 

Le traitement de données est normalement beaucoup plus rapide car les outils `bowtie2-build`, `bowtie` et `samtools` utilisent plusieurs coeurs simultanément.


@JULIEN :

- je voudrais maintenant utiliser l'option `--threads` de bowtie2 et de samtools qui permettent d'utiliser plusieurs threads. Je souhaite toujours n'utiliser qu'un seul noeud mais plusieurs coeurs de ce noeud. Ce que j'ai indiqué dans dans l'entête du script (`--cpus-per-task=8`) et au niveau de la commande `srun` qui lance bowtie2 et samtools te semble correct ?
