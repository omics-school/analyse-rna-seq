---
title: Calcul sur le cluster de l'IFB
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---


Dans cette activité, vous allez analyser les données RNA-seq de *O. tauri* sur le cluster NNCR de l'Institut Français de Bioinformatique (IFB). Ce cluster utilise un système d'exploitation Unix / Linux.

# Remarques préables

L'accès au cluster de l'IFB vous est fourni dans le cadre du DU Omique. C'est accès sera évoqué à l'issue de la formation, fin janvier 2020. Si vous souaitez continuer à utiliser ce cluster, n'hésitez pas à faire une demande en remplissant le formulaire [IFB core cluster - account request](https://www.france-bioinformatique.fr/fr/ifb-core-cluster-account-request) en précisant en quelques mots votre projet.

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

Quelles sont les versions de ces outils ? Si besoin, retournez voir le [Tutoriel de l'analyse RNA-seq](analyse_RNA-seq_O_tauri.md). Est-ce que sont les mêmes versions que sur le serveur du DU ?


## Préparation des données

Dans votre répertoire de travail (`/shared/projects/du_o_2019/login`), créez le répertoire `RNAseq`.

Copiez à l'intérieur les fichiers dont vous aurez besoin pour travailler :

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

FastQC va produire deux fichiers (un fichier avec l'extension `.html` et un autre avec l'extension `.zip`). Copiez le fichier `.html` sur votre machine locale avec le logiciel FileZilla ou la commande `scp`.


## Automatisation 1 

Téléchargez le script 3 avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script3.sh
```

Modifiez-le avec `nano` pour l'adapter à vos échantillons puis lancez-le avec la commande :
```
srun bash script3.sh
```

## Automatisation 2


Téléchargez le script 4 avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script4.sh
```


