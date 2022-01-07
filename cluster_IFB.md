---
title: Calcul sur le cluster de l'IFB
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---


Dans cette activité, vous allez analyser des données RNA-seq de *O. tauri* avec le cluster *National Network of Computational Resources* (NNCR) de l'[Institut Français de Bioinformatique]((https://www.france-bioinformatique.fr/) (IFB). Ce cluster utilise un système d'exploitation Linux.


## Remarques préables

L'accès au cluster de l'IFB vous est fourni dans le cadre pédagogique du DU Omiques. Cet accès sera révoqué à l'issue de la formation. 

Si, à l'issue de cette formation, vous souhaitez continuer à utiliser ce cluster pour votre projet de recherche, connectez-vous à l'[interface de gestion de votre compte IFB](https://my.cluster.france-bioinformatique.fr/manager2/project) puis cliquez sur le bouton *Request A New Project* et détaillez en quelques mots votre projet. Plusieurs utilisateurs peuvent être associées à un même projet et partager des données. Selon la quantité de ressources que vous demanderez, la création d'un projet pourra être associée à un coût. Au 06/01/2022, la grille tarifaire n'est pas encore connue.

Si vous avez besoin d'un logiciel spécifique sur le cluster. N'hésitez pas à le demander gentillement sur le site [Cluster Community Support](https://community.france-bioinformatique.fr/c/ifb-core-cluster/). Les administrateurs sont très sympas et en général très réactifs.


## 0. Connexion au cluster

Sous Windows, ouvrez un terminal Ubuntu. Si vous avez oublié comment faire, consultez le [tutoriel Unix](https://omics-school.github.io/unix-tutorial/tutoriel/README).

Connectez-vous en SSH au cluster avec les identifiants (login et mot de passe) que vous avez du recevoir par e-mail.

La syntaxe est de la forme :

```bash
$ ssh LOGIN@core.cluster.france-bioinformatique.fr
```

avec `LOGIN` votre identifiant sur le cluster. 

🔔 Rappels : 

- Ne tapez pas le `$` en début de ligne et faites attention aux majuscules et aux minuscules !
- Utilisez le copier / coller.
- Utilisez la complétion des noms de fichier et de répertoires avec la touche <kbd>Tab</kbd>.

Si c'est la première fois que vous vous connectez au cluster, répondez `yes` à la question 

```
Are you sure you want to continue connecting (yes/no)?
```

puis validez.

Entrez ensuite votre mot de passe en **aveugle**, c'est-à-dire sans qu'aucun caractère ne soit affiché à l'écran. C'est assez déstabilisant la première fois puis on s'habitue.

Si vous le souhaitez, une [vidéo](https://youtu.be/7OnrlZbVtEk) illustrant pas-à-pas la connexion en SSH est disponible.

🔔 **Attention** 🔔 Le cluster est protégé contre certaines attaques. Si vous entrez un mot de passe erronné plusieurs fois de suite, votre adresse IP va être bannie et vous ne pourrez plus vous connecter (temporairement) au serveur.

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

Utilisez la commande suivante pour compter le nombre de logiciels et de versions disponibles avec `module` :

```bash
$ module avail -l | wc -l
```

Vérifiez maintenant si la version 2.3.6 de `bowtie2` est disponible avec la commande :

```bash
$ module avail -l bowtie2
```

Si un jour vous avez besoin d'un logiciel dans une version spécifique, n'hésitez pas à le demander au [support communautaire](https://community.france-bioinformatique.fr/c/ifb-core-cluster/) du cluster.

Chargez ensuite les logiciels `fastqc`, `bowtie2`, `samtools` et `htseq` avec les commandes suivantes :

```bash
$ module load fastqc/0.11.9
$ module load bowtie2/2.3.5
$ module load samtools/1.9
$ module load htseq/0.11.3
```

**Remarque** : les logiciels chargés avec `module` ne sont disponibles que le temps de votre session sur le cluster. Si vous vous déconnectez puis vous vous reconnectez, il faudra charger à nouveau les logiciels dont vous aurez besoin.


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

Votre répertoire utilisateur sur le noeud de connexion est : `/shared/home/LOGIN` (avec `LOGIN` votre identifiant sur le cluster).

Ce répertoire ne doit pas contenir de donnée volumineuse car l'espace disponible est limité à 100 Go. Un espace de stockage plus conséquent a été créé pour vous dans le répertoire  `/shared/projects/form_2021_29/LOGIN` . Par la suite, cet espace sera appelé « répertoire de travail ».

De plus, le répertoire `/shared/projects/form_2021_29/data/rnaseq_tauri` contient les données dont vous aurez besoin pour ce projet. Vous n'avez accès à ce répertoire qu'en lecture seule, c'est-à-dire que vous pouvez seulement parcourir les répertoires et lire les fichiers de ce répertoire (pas de modification, d'ajout ou de suppression).

Essayez de vous souvenir de quels fichiers vous aurez besoin pour l'analyse des données RNA-seq de *O. tauri*. Trouvez-vous ces fichiers dans le répertoire `/shared/projects/form_2021_29/data/rnaseq_tauri` ou un de ses sous-répertoires ?

Vérifiez l'intégrité des fichiers `.fastq.gz` situés dans le répertoire `/shared/projects/form_2021_29/data/rnaseq_tauri/reads` avec les commandes suivantes :

```bash
$ cd /shared/projects/form_2021_29/data/rnaseq_tauri
```

*Rappel : N'entrez pas le symbole $ en début de ligne.*

puis :

```bash
$ srun -A form_2021_29 md5sum -c reads_md5sum.txt
```

N'oubliez pas `srun -A form_2021_29` en début de commande :

- L'instruction `srun` est spécifique au cluster. 
- L'option `-A form_2021_29` spécifie quel projet utiliser (facturer) pour cette commande. Un même utilisateur peut appartenir à plusieurs projets. Le nombre d'heures de calcul attribuées à un projet étant limité, il est important de savoir quel projet imputer pour telle ou telle commande. Pensez-y pour vos futurs projets.


Déplacez-vous maintenant dans votre répertoire de travail `/shared/projects/form_2021_29/LOGIN` (avec `LOGIN` votre identifiant sur le cluster). Un moyen simple d'y parvenir est d'exécuter la commande :

```bash
$ cd /shared/projects/form_2021_29/$USER
```

Créez ensuite le répertoire `rnaseq_tauri` et déplacez-vous à l'intérieur. Dorénavant vous ne travaillerez plus qu'à partir de ce répertoire.

La commande `pwd` devrait vous renvoyer quelque chose du type :

```bash
$ pwd
/shared/projects/form_2021_29/LOGIN/rnaseq_tauri`
```

avec `LOGIN` votre identifiant sur le cluster.

🛑 Si vous ne parvenez pas à être dans le bon répertoire, n'allez pas plus loin et appelez à l'aide 🆘.


## 3.1 Analyse d'un seul échantillon

En guise d'introduction, vous affichez tous les jobs en cours d'exécution sur le cluster avec la commande :

```bash
$ squeue -t RUNNING
```

Vous voyez que vous n'êtes pas seul ! Comptez maintenant le nombre de jobs en cours d'exécution en chaînant la commande précédente avec `wc -l` :

```bash
$ squeue -t RUNNING | wc -l
```


**Remarques importantes concernant l'indexation du génome de référence** : 

- L'indexation du génome de référence avec le logiciel `bowtie2` a déjà été effectué pour vous. Pour vous en convraincre, affichez le contenu du répertoire `/shared/projects/form_2021_29/data/rnaseq_tauri/genome` et vérifiez l'existence de fichiers avec l'extension `.bt2`, spécifiques des fichiers index créés par `bowtie2`.
- Cette indexation a été réalisée avec la commande `sbatch -A form_2021_29 /shared/projects/form_2021_29/data/rnaseq_tauri/build_genome_index.sh` que, bien sûr, vous n'exécuterez pas !

Depuis le cluster de l'IFB, vérifiez que vous êtes toujours dans votre répertoire `/shared/projects/form_2021_29/LOGIN/rnaseq_tauri`.

Téléchargez ensuite le script `script4.sh` avec la commande :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script4.sh
```

Pour en comprendre le fonctionnement, explorez le contenu de ce script avec la commande `less` :

```bash
$ less script4.sh
```

Utilisez les touches <kbd>↓</kbd> et <kbd>↑</kbd> pour naviguer dans le fichier. Pressez la touche <kbd>Q</kbd> pour quitter `less`.

Lancez ensuite ce script avec la commande :

```bash
$ sbatch -A form_2021_29 script4.sh
```

Vous devriez obtenir un message du type `Submitted batch job 20716345`. Ici, `20716345` correspond au numéro du job.

Notez bien le numéro de votre job.

Vérifiez que votre script est en train de tourner avec la commande :

```bash
$ squeue -u $USER
```

**Remarque** : Voici quelques statuts (colonne `ST`) de job intéressants :

- `CA` (*cancelled*) : le job a été annulé
- `F` (*failed*) : le job a planté
- `PD` (*pending*) : le job est en attente que des ressources soient disponibles
- `R` (*running*) : le job est en cours d'exécution


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

Si vous affichez le contenu de votre répertoire courant, vous devriez observer un fichier `slurm-JOBID.out` où `JOBID` est le numéro de votre job. Ce fichier contient la sortie, c'est-à-dire le *log* de votre script.

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

Faites aussi un peu de ménage en supprimant les fichiers créés précédemment :

```bash
$ rm -rf map/ reads_qc/ count/ slurm*.out
```

## 3.2 Analyse plus rapide d'un échantillon

L'objectif est maintenant « d'aller plus vite » en attribuant plusieurs coeurs pour l'étape d'alignement des reads sur le génome avec `bowtie2`.

Toujours depuis le cluster de l'IFB, dans le répertoire `rnaseq_tauri` de votre répertoire de travail, téléchargez le script `script5.sh` avec la commande :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script5.sh
```

Identifiez les différences avec le script précédent avec la commande `diff` :

```bash
$ diff script4.sh script5.sh
```

Les lignes qui débutent par `<` viennent de `script4.sh` et celles qui débutent par `>` viennent de `script5.sh`.

La différence majeure avec `script4.sh` réside dans l'utilisation de plusieurs coeurs pour la commande `bowtie2` avec l'option `--threads="${SLURM_CPUS_PER_TASK}"`. L'utilisation de plusieurs coeurs est permise par la déclaration `#SBATCH --cpus-per-task=8` au tout début de `script5.sh`.

**Remarque** : nous aurions également pu attribuer plusieurs coeurs pour les commandes `samtools view` et `samtools sort`, mais nos tests ont montré qu'il n'y avait pas, pour ce cas précis, de gain significatif en terme de temps de calcul. Pour information, les lignes de commande à utiliser seraient :

```bash
srun samtools view --threads="${SLURM_CPUS_PER_TASK}" -b "map/bowtie-${sample}.sam" -o "map/bowtie-${sample}.bam"
srun samtools sort --threads="${SLURM_CPUS_PER_TASK}" "map/bowtie-${sample}.bam" -o "map/bowtie-${sample}.sorted.bam"
```

Notez que tous les logiciels ne peuvent utiliser plusieurs coeurs. Consultez toujours la documentation de l'outil considéré.

Lancez maintenant le script d'analyse `script5.sh` :

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

- La commande `watch` est utilisée ici pour « surveiller » en quasi-temps réel le résultat de la commande `sacct`.
- L'affichage est rafraichi toutes les 2 secondes.
- Vous pouvez également afficher la mémoire vive maximale consommée à chaque étape avec la commande `sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,MaxRSS,NodeList -j JOBID`

Votre job devrait prendre une dizaine de minutes pour s'exécuter. Laissez le cluster travailler et profitez-en pour vous préparer un thé ou un café bien mérité.

Quand les statuts (colonne `State`) du job et de tous ses « *job steps* » sont à `COMPLETED`, quittez la commande `watch` en appuyant sur la combinaison de touches <kbd>Ctrl</kbd> + <kbd>C</kbd>.

Vérifiez avec la commande `tree` que vous obtenez une arborescence équivalente à celle ci-dessous :

```
$ tree
.
├── count
│   └── count-SRR2960338.txt
├── map
│   ├── bowtie-SRR2960338.sorted.bam
│   └── bowtie-SRR2960338.sorted.bam.bai
├── reads_qc
│   ├── SRR2960338_fastqc.html
│   └── SRR2960338_fastqc.zip
├── script4.sh
├── script5.sh
└── slurm-20716384.out
```

Vérifiez que la somme de contrôle du fichier `count/count-SRR2960338.txt` est bien `36fc86a522ee152c89fd77430e9b56a5`.

Faites maintenant un peu de ménage en supprimant les fichiers créés précédemment avec la commande :

```bash
$ rm -rf map/ reads_qc/ count/ slurm*.out
```


## 3.3 Analyse de plusieurs échantillons

Toujours depuis le cluster de l'IFB, dans le répertoire `rnaseq_tauri` de votre répertoire de travail, téléchargez le script `script6.sh` avec la commande :

```bash
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script6.sh
```

Nous pourrions analyser d'un seul coup les 47 échantillons (fichiers `.fastq.gz`) mais pour ne pas consommer trop de ressources sur le cluster, nous allons limiter notre analyse à 4 échantillons seulement. Si vous le souhaitez vous pourrez modifier ce script pour analyser les 47 échantillons 💪.

Lancez votre analyse avec la commande :

```bash
$ sbatch -A form_2021_29 script6.sh
```

Notez bien le numéro de job renvoyé.

Suivez en temps réel l'exécution de votre job avec la commande :

```bash
$ watch sacct --format=JobID,JobName,State,Start,Elapsed,CPUTime,NodeList -j JOBID
```
avec `JOBID` le numéro de votre job.

Remarquez que la ligne indiquant `script6.sh` pour « *JobName* » est présente 4 fois. Cela indique que le script `script6.sh` est exécutée 4 fois, en parallèle.

Patientez une dizaine de minutes que tous les jobs et *job steps* soient terminés. 

Quand les status (colonne `State`) de tous les jobs et *job steps* sont à `COMPLETED`, quittez la commande `watch` en appuyant sur la combinaison de touches <kbd>Ctrl</kbd> + <kbd>C</kbd>.

**Remarques**: 

- Notez que l'exécution de `script6.sh` aura pris environ le même temps que celle de `script5.sh`. C'est toute la puissance du calcul distribué 🚀. 
- Vous comprenez qu'il est possible d'analyser 4, 10 ou 47 échantillons dans un temps raisonnable. Cependant, si vous lancez une telle analyse sur un grand nombre d'échantillons, il est possible que tous vos jobs ne partent pas en même temps et que certains aient le statut *PENDING*, le temps que les ressources nécessairent se libèrent sur le cluster.
- L'analyse de la totalité des 47 échantillons (fichiers .fastq.gz) génère environ 18 Go de données.


Une dernière fois, vérifiez que tous vos fichiers sont présents dans les bons répertoires :

```bash
$ tree
.
├── count
│   ├── count-SRR2960338.txt
│   ├── count-SRR2960341.txt
│   ├── count-SRR2960343.txt
│   └── count-SRR2960356.txt
├── map
│   ├── bowtie-SRR2960338.sorted.bam
│   ├── bowtie-SRR2960338.sorted.bam.bai
│   ├── bowtie-SRR2960341.sorted.bam
│   ├── bowtie-SRR2960341.sorted.bam.bai
│   ├── bowtie-SRR2960343.sorted.bam
│   ├── bowtie-SRR2960343.sorted.bam.bai
│   ├── bowtie-SRR2960356.sorted.bam
│   └── bowtie-SRR2960356.sorted.bam.bai
├── reads_qc
│   ├── SRR2960338_fastqc.html
│   ├── SRR2960338_fastqc.zip
│   ├── SRR2960341_fastqc.html
│   ├── SRR2960341_fastqc.zip
│   ├── SRR2960343_fastqc.html
│   ├── SRR2960343_fastqc.zip
│   ├── SRR2960356_fastqc.html
│   └── SRR2960356_fastqc.zip
├── script4.sh
├── script5.sh
├── script6.sh
├── slurm-20716400_0.out
├── slurm-20716400_1.out
├── slurm-20716400_2.out
└── slurm-20716400_3.out
```

Comme vous avez lancé 4 sous-jobs indépendants, SLURM a également créé 4 fichiers de sortie distincts.

## 4. L'heure de faire les comptes

Expérimentez la commande `sreport` pour avoir une idée du temps de calcul consommé par tous vos jobs :

```bash
$ sreport -t hour Cluster UserUtilizationByAccount Start=2022-01-01 End=$(date --iso-8601)T23:59:59 Users=$USER
```

La colonne `Used` indique le nombre d'heures de temps CPU consommées. Cette valeur est utile pour estimer le « coût CPU » d'un projet.

Voici un exemple de rapport produit par `sreport` :

```bash
$ sreport -t hour Cluster UserUtilizationByAccount Start=2022-01-01 End=$(date --iso-8601)T23:59:59 Users=$USER
--------------------------------------------------------------------------------
Cluster/User/Account Utilization 2022-01-01T00:00:00 - 2022-01-06T17:59:59 (496800 secs)
Usage reported in CPU Hours
--------------------------------------------------------------------------------
  Cluster     Login     Proper Name         Account     Used   Energy 
--------- --------- --------------- --------------- -------- -------- 
     core  ppoulain  Pierre Poulain    form_2021_29       93        0 
     core  ppoulain  Pierre Poulain          gonseq        5        0 
```

Ainsi, l'utilisateur `ppoulain` a déjà consommé 93 heures de temps CPU sur le projet `form_2021_29`.

Attention, `sreport` ne prend pas en compte les heures immédiatement consommées. Il lui faut un peu de temps pour consolider les données.


## 5. Quelques conseils supplémentaires

L'analyse RNA-seq présentée ici tourne en 10-15', c'est très rapide car le génome d'*O. tauri* est relativement petit. Les temps d'analyse seront plus longs avec des génomes plus gros.

Procédez toujours par itérations successives. Testez votre script d'analyse RNA-seq pour 1 échantillon, puis 2 ou 3 puis la totalité.

Quand vous lancez un job qui sera potentiellement long, n'hésitez pas à ajouter les directives ci-dessous au début de votre script :

```
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=votre-adresse-mail@email.fr
```

Vous recevrez alors automatiquement un e-mail lorsque le job se termine ou si celui-ci plante.

Sur le cluster IFB :

- Un utilisateur ne peut utiliser plus de 300 coeurs en même temps.
- Un job dure, par défaut, au maximum 24 h. Une queue plus longue (appelée `long`) est disponible pour des jobs qui durent jusqu'à 30 jours et est utilisable via la directive `#SBATCH --partition=long` en début de script. D'autres queues plus spécifiques (plus de temps, beaucoup de mémoire vive, GPU) sont disponibles sur demande.

Prenez le temps d'explorer la [documentation très complète](https://ifb-elixirfr.gitlab.io/cluster/doc/) sur le cluster. Vous y trouverez notamment un [tutoriel](https://ifb-elixirfr.gitlab.io/cluster/doc/tutorials/analysis_slurm/) sur une autre analyse RNA-seq.


## 6. Sur place ou à emporter ?

Nous allons maintenant utiliser deux stratégies pour récupérer les résultats de l'analyse RNA-seq et copier les données depuis le cluster vers notre machine locale.

Une [vidéo](https://youtu.be/1UTPzK-zeA8) illustrant ces deux approches est également disponible.

### 6.1 scp

⚠️ Pour récupérer des fichiers sur le cluster en ligne de commande, vous devez lancer la commande `scp` depuis un shell Unix sur votre machine locale. ⚠️

Depuis un shell Unix sur votre machine locale, déplacez-vous dans le répertoire `/mnt/c/Users/omics` et créez le répertoire `rnaseq_tauri_cluster`. 

Déplacez-vous dans ce nouveau répertoire.

Utilisez la commande `pwd` pour vérifier que vous êtes bien dans le répertoire `/mnt/c/Users/omics/rnaseq_tauri_cluster`. 

Lancez ensuite la commande suivante pour récupérer les fichiers de comptage :

```bash
$ scp LOGIN@core.cluster.france-bioinformatique.fr:/shared/projects/form_2021_29/LOGIN/rnaseq_tauri/count/count*.txt .
```

où `LOGIN` est votre identifiant sur le cluster (qui apparait deux fois dans la ligne de commande ci-dessus). Faites bien attention à garder le `.` à la fin de la ligne de commande.

Comme d'habitude, entrez votre mot de passe du cluster en aveugle.

Une fois que vous avez récupéré les résultats de comptage, vérifiez que la somme de contrôle MD5 du fichier `count-SRR2960338.txt` est bien la même que précédemment (`36fc86a522ee152c89fd77430e9b56a5`).

Pour récupérer directement le répertoire `count` sur le cluster, vous auriez pu utiliser la commande :

```bash
$ scp -r LOGIN@core.cluster.france-bioinformatique.fr:/shared/projects/form_2021_29/LOGIN/rnaseq_tauri/count .
```

Notez l'option `-r` qui indique qu'on transfère un répertoire.

### 6.2 FileZilla

Lancez le logiciel FileZilla ([comme ceci](img/filezilla.png)). Puis entrez les informations suivantes :

- Hôte : `sftp://core.cluster.france-bioinformatique.fr`
- Identifiant : votre login sur le cluster
- Mot de passe : votre mot de passe sur le cluster

Cliquez ensuite sur le bouton *Connexion rapide*. Cliquez sur *OK* dans la fenêtre *Clé de l'hôte inconnue*

Une fois connecté :

- Dans le champ texte à côté de *Site local* (à gauche de la fenêtre), entrez le chemin `C:\Users\omics\rnaseq_tauri_cluster`.
- Dans le champ texte à coté de *Site distant* (à droite de la fenêtre), entrez le chemin `/shared/projects/form_2021_29/` voire directement votre répertoire de travail `/shared/projects/form_2021_29/LOGIN/` (avec `LOGIN` votre identifiant sur le cluster).

Essayez de transférer des fichiers dans un sens puis dans l'autre. Double-cliquez sur les fichiers pour lancer les transferts.



