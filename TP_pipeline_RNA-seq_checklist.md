---
title: Check list analyse de données RNA Seq
author: Pierre Poulain
license: Creative Commons Attribution (CC-BY)
---

# Analyse de données RNA-seq sous Linux : check list

## Préparation de l'environnement de travail

- [ ] Conda est disponible dans mon shell.  
    La commande ```conda --version``` ne renvoie pas de message d'erreur.  
    La version de conda est : ___________________________
- [ ] L'environnement conda pour l'analyse RNA-seq est activé.  
    Mon invite de commande ressemble à quelque chose du type :
    ```
    (rnaseq) ppoulain@candihub:~$
    ```
- [ ] La version de `fastqc` est : ___________________________
- [ ] La version de `bowtie2` est : ___________________________
- [ ] La version de `samtools` est : ___________________________
- [ ] La version de `htseq-count` est : ___________________________


## Préparation des données

- [ ] Le répertoire `RNAseq` est créé dans mon répertoire personnel.
- [ ] Les fichiers contenant les reads (`.fastq.gz`) sont copiés dans le répertoire `RNAseq`.  
    Identifiants des échantillons (le `XX` de `140317_SN365_A_L001_HCA-XX_R1.fastq.gz`)  
    : __________________________________________________
- [ ] Le fichier contenant la séquence du génome de référence (`.fna`) est copié dans le répertoire `RNAseq`.
- [ ] Le fichier contenant les annotations du génome de référence (`_DUO2.gff`) est copié dans le répertoire `RNAseq`.


## Analyse manuelle

Identifiant d'un échantillon contenant les *reads* : ___

### Contrôle qualité

- [ ] Lancement de FastQC.
- [ ] Copie du fichier `.html` généré par FastQC sur votre machine (via FileZilla ou scp).
- [ ] Visualisation de l'analyse.


### Indexation du génome de référence

- [ ] Lancement de Bowtie (commande `bowtie2-build`).
- [ ] Taille occupée par les fichiers d'index crées par Bowtie : ____


### Alignements des *reads* sur le génome de référence

- [ ] Lancement de Bowtie (commande `bowtie2`).  
    Nombre total de *reads* dans le fichier `.fastq.gz` : ____________  
    Nombre de *reads* non alignés : ____________  
    Nombre de *reads* alignés (au moins 1 fois) : ____________  
    Taux d'alignement : ____________  
- [ ] Taille du fichier d'alignement créé par Bowtie : ___________

Pourquoi le fichier d'alignement créé par Bowtie est beaucoup plus gros que le fichier contenant les *reads* ?


### Conversion des *reads* alignés en binaire, tri et indexation

- [ ] Lancement de SAMtools pour convertir le fichier d'alignement des *reads* en binaire.
    Taille du fichier créé par samtools : _______________
- [ ] Lancement de samtools pour trier les *reads* alignés.
    Taille du fichier créé par samtools : _______________
- [ ] Lancement de samtools pour indexer les *reads* alignés.
    Nom du fichier d'index : ___________________________________


### Visualisation des *reads* alignés avec IGV

Les fichiers suivants sont copiés sur la machine locale :
- [ ] Génome de référence (`.fna`)
- [ ] Annotations du génome de référence (`_DUO2.gff`)
- [ ] bam trié (`bowtie.sorted.bam`)
- [ ] index du bam trié (`bowtie.sorted.bam.bai`)

IGV :
- [ ] Lancement d'IGV et chargement des différents fichiers
- [ ] Visualisation du gène `ostta18g01980` dans IGV.


### Comptage des *reads* alignés sur les gènes de *O. tauri*

- [ ] Lancement de HTSeq
    Nombre d'annotations contenues dans le fichier `.gff` : _______________
- [ ] Recherche du gène `ostta18g01980` dans le fichier de comptage
    Nombre de *reads* alignés sur ce gène : _______________

Bravo !


## Automatisation de l'analyse : niveau 1



## Automatisation de l'analyse : niveau 2



## Automatisation de l'analyse : niveau 3
