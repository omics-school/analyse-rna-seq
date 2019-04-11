---
title: Check list analyse de données RNA-seq O. tauri avec Unix
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

# Analyse de données RNA-seq avec Unix : *check list*

## Étape 1 : préparation de l'environnement de travail

### Chargement de l'environnement conda

- [ ] Conda est disponible dans mon shell.  
    La commande ```conda --version``` ne renvoie pas de message d'erreur.  
    La version de conda est : ___________________________
- [ ] L'environnement conda pour l'analyse RNA-seq est activé.  
    Mon invite de commande ressemble à quelque chose du type :
    ```
    (rnaseq) ppoulain@candihub:~$
    ```

### Vérification des logiciels disponibles

- [ ] La version de `fastqc` est : ___________________________
- [ ] La version de `bowtie2` est : ___________________________
- [ ] La version de `samtools` est : ___________________________
- [ ] La version de `htseq-count` est : ___________________________


## Étape 2 : préparation des données

- [ ] Le répertoire `RNAseq` est créé dans mon répertoire personnel.
- [ ] Les fichiers contenant les *reads* (`.fastq.gz`) sont copiés dans le répertoire `RNAseq`.  
    Identifiants des échantillons (le `XX` de `140317_SN365_A_L001_HCA-XX_R1.fastq.gz`)  
    ._______ _______ _______ _______
- [ ] Le fichier contenant la séquence du génome de référence (`.fna`) est copié dans le répertoire `RNAseq`.
- [ ] Le fichier contenant les annotations du génome de référence (`.gff`) est copié dans le répertoire `RNAseq`.
- [ ] Les fichiers contenant les *reads* (`.fastq.gz`) ont été renommés correctement sur la forme `HCA-XX_R1.fastq.gz`, avec `XX` le numéro de l'échantillon.
- [ ] Volume total des données copiées : ___________


## Étape 3 : analyse manuelle

Identifiant de l'échantillon contenant les *reads* que vous allez analyser : ___

### Contrôle qualité

- [ ] Lancement de FastQC.
- [ ] Copie du fichier `.html` généré par FastQC sur ma machine (via FileZilla ou scp).
- [ ] Visualisation de l'analyse.


### Indexation du génome de référence

- [ ] Lancement de Bowtie (commande `bowtie2-build`).
- [ ] Taille occupée par les fichiers d'index crées par Bowtie : ______


### Alignements des *reads* sur le génome de référence

- [ ] Lancement de Bowtie (commande `bowtie2`).  
    Nombre total de *reads* dans le fichier `.fastq.gz` : ____________  
    Nombre de *reads* non alignés : ____________  
    Nombre de *reads* alignés (au moins 1 fois) : ____________  
    Taux d'alignement : ____________  
- [ ] Taille du fichier d'alignement créé par Bowtie : ___________

Pourquoi le fichier d'alignement créé par Bowtie est beaucoup plus gros que le fichier contenant les *reads* ?  
.______________________________________________________________
.______________________________________________________________
.______________________________________________________________


### Conversion des *reads* alignés en binaire, tri et indexation

- [ ] Lancement de SAMtools pour convertir le fichier d'alignement des *reads* en binaire.  
    Taille du fichier créé par SAMtools : _______________
- [ ] Lancement de samtools pour trier les *reads* alignés.  
    Taille du fichier créé par SAMtools : _______________
- [ ] Lancement de samtools pour indexer les *reads* alignés.  
    Nom du fichier d'index : ___________________________________


### Visualisation des *reads* alignés avec IGV

Les fichiers suivants sont copiés sur la machine locale :
- [ ] génome de référence (`.fna`)
- [ ] annotations du génome de référence (`_DUO2.gff`)
- [ ] bam trié (`bowtie.sorted.bam`)
- [ ] index du bam trié (`bowtie.sorted.bam.bai`)

IGV :
- [ ] Lancement d'IGV et chargement des différents fichiers
- [ ] Visualisation du gène `ostta18g01980` dans IGV.


### Comptage des *reads* alignés sur les gènes de *O. tauri*

- [ ] Lancement de HTSeq  
    Nombre d'annotations contenues dans le fichier `.gff` : _______________
- [ ] Recherche du gène `ostta18g01980` dans le fichier de comptage.  
    Nombre de *reads* alignés sur ce gène : _______________

Bravo !


## Étape 4 : automatisation de l'analyse : niveau 1

- [ ] Ouverture du script 1
- [ ] Les trois variables dans ce script sont :
    ._______________ _______________ _______________
- [ ] Téléchargement du script 1
- [ ] Modification de la variable qui contient le numéro d'échantillon
- [ ] Lancement du script 1
- [ ] Café ! ☕ 🍪 ☕ 🍪


## Étape 5 : automatisation de l'analyse : niveau 2

- [ ] Une évolution possible du premier script :
    ._____________________________________________
- [ ] Téléchargement du script 2
- [ ] Modification du script 2 pour mon échantillon.
- [ ] Lancement du script 2.


## Étape 6 : Automatisation de l'analyse : niveau 3 (ninja)

- [ ] Activité et exercices *boucle* de Software Carpentry
- [ ] Téléchargement du script 3
- [ ] La variable à modifier avec mes numéros d'échantillon est :
    ._________________
- [ ] Lancement du script 3.
