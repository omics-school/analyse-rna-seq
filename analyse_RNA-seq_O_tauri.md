---
title: Analyse des données RNA-seq O. tauri avec Unix
author: Pierre Poulain
license: Creative Commons Attribution-ShareAlike (CC BY-SA 4.0)
---

# Analyse des données RNA-seq de *O. tauri* avec Unix

Dans cette activité, vous allez analyser les données RNA-seq de *O. tauri* dans un environnement Unix/Linux.

Pour cela, vous allez beaucoup utiliser la ligne de commande, connecté en SSH sur le serveur du DU. Vous copierez également des fichiers depuis le serveur du DU vers votre machine locale avec le logiciel FileZilla ou la commande `scp` (nous y reviendrons) .

Voici une vue d'ensemble des étapes pour analyser les données de séquençage haut débit :

![](pipeline_RNA_seq_O_tauri.svg)


## Étape 1 : préparation de l'environnement de travail

Pour analyser les données de séquençage haut débit de *O. tauri.*, nous avons besoin des logiciels suivants :

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) pour le contrôle qualité des données de séquençage.
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) pour l'indexation du génome de référence puis l'alignement des *reads* sur le génome de référence.
- [SAMtools](http://samtools.sourceforge.net/) pour la manipulation des fichiers d'alignements (conversion en binaire, tri et indexation)
- [HTSeq](https://htseq.readthedocs.io/en/latest/) pour le comptage du nombre de *reads* alignés sur chaque gène.


### Configuration de conda

La distribution Miniconda a été installée sur le serveur du DU. Pour que vous puissiez avoir accès à cet outil,
vous devez configurer votre *shell* Linux sur le serveur du DU. Les étapes à suivre sont :

1. Connectez-vous en SSH au serveur du DU.
1. Éditez le fichier `.bashrc` dans votre répertoire personnel. Par exemple avec l'éditeur de texte nano :
    ```
    $ nano ~/.bashrc
    ```
1. Déplacez-vous à la fin du fichier et ajoutez la ligne ci-dessous :
    ```
    source /data/omics-school/share/miniconda/etc/profile.d/conda.sh
    ```
    Enregistrez le fichier avec la combinaison de touches <kbd>Ctrl</kbd> + <kbd>O</kbd> puis validez par <kbd>Entrée</kbd>.
    Puis quittez nano avec la combinaison de touches <kbd>Ctrl</kbd> + <kbd>X</kbd>.  
    Remarque 1 : la ligne de commande à ajouter est assez longue. Pour éviter les erreurs, utilisez le copier (<kbd>Ctrl</kbd> + <kbd>C</kbd>) / coller (clic droit) dans nano.  
    Remarque 2 : il est possible que votre fichier `.bashrc` soit vide, ce n'est pas un problème.
1. Vérifiez que conda est maintenant disponible en vous déconnectant du serveur, en vous reconnectant puis en tapant la commande suivante :
    ```
    $ conda --version
    ```
1. Bravo ! 🎉 Vous avez correctement configuré conda.

Les manipulations ci-dessus vous ont permis de rendre disponible conda dans votre *shell* Linux sur le serveur du DU. Elles ne sont à faire qu'une seule fois.

Une documentation expliquant l'installation de miniconda et la configuration de conda est disponible [ici](conda.md).


### Chargement de l'environnement conda

Conda est un gestionnaire de paquets qui permet d'installer de nombreux logiciels utilisés en bioinformatique. Il permet aussi de créer des environnements *virtuels* dans lequel ces logiciels sont installés. L'intérêt des environnements virtuels est de pouvoir installer sur la même machine, plusieurs version d'un même logiciel, chaque version étant installée dans un environnement virtuel indépendant.

Nous allons maintenant voir comment charger un environnement virtuel créé avec conda.

1. Connectez-vous en SSH au serveur du DU.
1. Vérifiez que conda est bien disponible avec la commande
    ```
    $ conda --version
    ```
1. Nous avons préparé un environnement virtuel conda spécialement pour l'analyse RNA-seq. Cet environement s'appelle `rnaseq`. Chargez cet environnement :
    ```
    $ conda activate rnaseq
    ```
1. Votre invite de commande devrait être modifiée et ressembler à :
    ```
    (rnaseq) ppoulain@candihub:~$
    ```
    La mention `(rnaseq)` au début de l'invite de commande indique que vous êtes dans l'environnement virtuel `rnaseq`.

Remarque : pour quitter un environnement virtuel, il faut utiliser la commande
```
$ conda deactivate
```
mais vous n'aurez pas besoin de l'utiliser pour cette activité.

Pour la suite, nous supposerons que :
1. Vous êtes connecté en SSH au serveur du DU.
1. Vous avez activé l'environnement conda `rnaseq`.


### Vérification des logiciels disponibles

En bioinformatique, il est essentiel de vérifier et noter les versions des logiciels que vous utilisez. Vous reporterez les noms et les versions des logiciels que vous utiliserez dans la section *Materials & Methods* de vos articles.

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
Certains programmes peuvent en effet renvoyer beaucoup d'informations.


### Comparaison avec les logiciels utilisés dans Galaxy (si vous avez du temps)

Connectez-vous maintenant à votre compte sur Galaxy. Essayez de retrouver les versions des logiciels que vous utilisés (FastQC, Bowtie2, SAMtools, HTSeq).

Pour ce faire, dans votre *History*, cliquez sur le nom d'un résultat d'analyse, puis cliquez sur le petit i entouré (:information_source:) et lisez les informations de la section *Job Dependencies*.

Comparez les versions des logiciels disponibles dans Galaxy et installés sur le serveur du DU.


## Étape 2 : préparation des données

Sur le serveur du DU, les données brutes dont vous aurez besoin sont dans le répertoire `/data/omics-school/share/tauri_2019/`. Voici un aperçu de l'organisation et du contenu de ce répertoire :
```
$ tree -h /data/omics-school/share/tauri_2019/ 
/data/omics-school/share/tauri_2019/
├── [462K]  GCF_000214015.3_version_140606_genomic_DUO2.gff
├── [ 13M]  GCF_000214015.3_version_140606_genomic.fna
└── [4.0K]  reads
    ├── [507M]  140317_SN365_A_L001_HCA-10_R1.fastq.gz
    ├── [510M]  140317_SN365_A_L001_HCA-11_R1.fastq.gz
    ├── [374M]  140317_SN365_A_L001_HCA-12_R1.fastq.gz
    ├── [399M]  140317_SN365_A_L001_HCA-13_R1.fastq.gz
    ├── [375M]  140317_SN365_A_L001_HCA-14_R1.fastq.gz
    ├── [441M]  140317_SN365_A_L001_HCA-15_R1.fastq.gz
    ├── [440M]  140317_SN365_A_L001_HCA-16_R1.fastq.gz
    ├── [587M]  140317_SN365_A_L001_HCA-17_R1.fastq.gz
    ├── [531M]  140317_SN365_A_L001_HCA-18_R1.fastq.gz
    ├── [944M]  140317_SN365_A_L001_HCA-19_R1.fastq.gz
    ├── [459M]  140317_SN365_A_L001_HCA-20_R1.fastq.gz
    ├── [907M]  140317_SN365_A_L001_HCA-21_R1.fastq.gz
    ├── [393M]  140317_SN365_A_L001_HCA-22_R1.fastq.gz
    ├── [429M]  140317_SN365_A_L001_HCA-23_R1.fastq.gz
    ├── [930M]  140317_SN365_A_L001_HCA-24_R1.fastq.gz
...
```

Dans votre répertoire personnel, créez le répertoire `RNAseq`. Faites attention aux minuscules et aux majuscules !

Dans ce répertoire `RNAseq`, copiez :

- Les 2 ou 3 fichiers contenant les *reads* (`.fastq.gz`) qui vous devez analyser. Tous les fichiers sont dans le répertoire  `/data/omics-school/share/tauri_2019/reads`
- Le génome de référence de *O. tauri* :  
    `/data/omics-school/share/tauri_2019/GCF_000214015.3_version_140606_genomic.fna`
- Les annotations du génome de référence :  
    `/data/omics-school/share/tauri_2019/GCF_000214015.3_version_140606_genomic_DUO2.gff`

Remarque : le génome de référence de *Ostreococcus tauri* et ses annotations sont disponibles sur la [page dédiée](https://www.ncbi.nlm.nih.gov/genome/373?genome_assembly_id=352933) sur le site du NCBI :
- [génome de référence](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.fna.gz)
- [annotations](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_genomic.gff.gz). Nous avons légèrement modifié le fichier d'annotations pour ne prendre en compte que les gènes et alléger la visualisation dans IGV.


**⚠ Étape essentielle pour la suite ⚠**

Renommez les fichiers contenant vos *reads* (`.fastq.gz`) de la façon suivante :
```
HCA-XX_R1.fastq.gz
```
avec `XX` le numéro de votre échantillon. Cela revient à supprimer la première partie du nom du fichier.

Par exemple, pour les fichiers correspondants aux échantillons 10 et 41 :
```
$ mv 140317_SN365_A_L001_HCA-10_R1.fastq.gz HCA-10_R1.fastq.gz
$ mv 140317_SN365_A_L002_HCA-41_R1.fastq.gz HCA-41_R1.fastq.gz
```

Cette opération simplifiera considérablement l'automatisation des analyses par la suite.

Déterminez la taille de vos données avec la commande 
```
$ du -ch *
```

Explications : la commande `du` affiche la taille occupée par des fichiers. L'option `-h` affiche la taille en ko, Mo, Go... L'option `-c` calcule la taille totale occupée par tous les fichiers.


## Étape 3 : analyse manuelle

Pour cette première analyse, choisissez un **seul échantillon** contenant des *reads*.


### Contrôle qualité

Lancez FastQC avec la commande :

```
$ fastqc nom-fichier-fastq.gz
```
où `nom-fichier-fastq.gz` est le fichier contenant l'échantillon que vous avez choisi.

FastQC va produire deux fichiers (un fichier avec l'extension `.html` et un autre avec l'extension `.zip`). Copiez le fichier `.html` sur votre machine locale avec le logiciel FileZilla ou la commande `scp`.


#### `scp`

Pour `scp`, vous devez être dans un *shell* sur votre **machine locale** et taper la commande 

```
$ scp login@omics-school.net:~/RNAseq/HCA-numéro_R1_fastqc.html ./
```
où `login` est votre identifiant sur le serveur du DU et `numéro` est le numéro de l'échantillon que vous avez analysé.

Entrez votre mot de passe lorsqu'on vous le demande.


#### FileZilla 

FileZilla est un logiciel open source (donc libre et gratuit) qui permet de transférer graphiquement des fichiers entre un serveur et une machine locale.

Installez tout d'abord FileZilla en le [téléchargeant](https://filezilla-project.org/).

Lancez-le puis remplissez les champs suivants : 

- Hôte : `sftp://omics-school.net`
- Identifiant : `<votre-login-sur-le-serveur>`
- Mot de passe : `<votre-mot-de-passe-sur-le-serveur>`

Puis cliquez sur le bouton *Connexion rapide*. 

Vous devriez obtenir sur le panneau de gauche l'arborescence de vos fichiers locaux et sur le panneau de droite l'arboresence de vos fichiers sur le serveur. Vous pouvez transférer les fichiers de l'un vers l'autre par glisser/déposer.


#### Analyse des résultats de FastQC

Une fois les bons fichiers transférés sur votre machine locale :

- Ouvrez ce fichier dans un navigateur internet (Firefox par exemple).
- Analysez le rapport de FastQC.


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


### Alignements des *reads* sur le génome de référence

Lancez l'alignement :
```
$ bowtie2 -x O_tauri -U nom-du-fichier.fastq.gz -S bowtie.sam
```

Ici :
- `O_tauri` désigne les fichiers index du génome de référence,
- `nom-fichier-fastq.gz` est le fichier contenant l'échantillon que vous avez choisi
- et `bowtie.sam` est le fichier qui va contenir l'alignement produit par Bowtie2.

Cette étape est la plus longue et peut prendre une dizaine de minutes. Bowtie n'affiche rien à l'écran lorsqu'il fonctionne. Soyez patient.

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

Pour visualiser l'alignement des *reads* sur le génome de référence avec IGV, copiez, avec FileZilla ou la commande `scp`, sur votre machine locale les fichiers :
- génome de référence (fichier `.fna`) ;
- annotations du génome de référence (fichier `_DUO2.gff`) ;
- bam trié (`bowtie.sorted.bam`) ;
- index du bam trié (`bowtie.sorted.bam.bai`).

Lancez IGV et visualisez l'alignement des *reads* sur le génome de référence. Si vous avez oublié comme faire, visionnez la vidéo sur ce sujet qui vous a été proposée précédemment.

Visualisez particulièrement le gène `ostta18g01980`.


### Comptage des *reads* alignés sur les gènes de *O. tauri*

Le comptage des *reads* alignés sur les gènes se fait avec HTSeq.

De retour sur le serveur, lancez la commande :
```
$ htseq-count --stranded=no --type='gene' --idattr='ID' --order=name --format=bam bowtie.sorted.bam GCF_000214015.3_version_140606_genomic_DUO2.gff > count.txt
```

HTSeq renvoie le nombre d'annotations trouvées dans le fichier `.gff` puis affiche une progression de l'analyse. Les options du programme `htseq-count` sont décrites dans la [documentation](http://gensoft.pasteur.fr/docs/HTSeq/0.9.1/count.html).

Déterminez le nombre de *reads* alignés sur le gène `ostta18g01980`. Pour cela, vous pouvez lancer la commande
```
$ grep ostta18g01980 count.txt
```
ou alors ouvrir le fichier `count.txt` avec la commande `less` puis chercher le gène `ostta18g01980` en tapant `/ostta18g01980` puis la touche `Entrée`.


### Trucs et astuces

Certaines étapes d'analyse (notamment l'alignement des *reads* sur le génome de référence et le comptage des *reads*) vont prendre du temps et consommer des ressources informatiques.

Si vous fermez votre terminal alors que vous avez une tache d'analyse en cours, celle-ci sera arrêtée. C'est dommage 😭

Pour lancer une analyse en tâche de fond et pouvoir vous déconnecter, utilisez la syntaxe :
```
$ nohup votre-commande-avec-ses-paramètres &
```
Attention, tout ce qui s'affiche normalement à l'écran sera écrit dans le fichier `nohup.out`. Vous obtiendrez d'ailleurs un message qui vous le confirme :
```
nohup: ignoring input and appending output to 'nohup.out'
```

Pour suivre l'avancée de votre analyse, lancez la commande
```
$ top
```
ou mieux, si le programme `htop` est installé sur le serveur :
```
$ htop
```

Enfin, voici quelques commandes utiles pour explorer les caractéristiques du serveur :
- `grep -c /proc/cpuinfo` pour connaître le nombre de coeurs de la machine
- `free -h` pour connaître la taille de la mémoire vive (total et disponible)
- `df -h` pour connaître la volumétrie de stockage (totale et diponible)
- `who` pour savoir qui est connecté sur le serveur


## Étape 4 : automatisation de l'analyse : niveau 1

Tout cela est très bien mais les fichiers que vous avez générés (`bowtie.bam`, `bowtie.sorted.bam`, `count.txt`...) ne sont pas très informatifs sur l'échantillon dont ils proviennent.

Par ailleurs, entrer toutes ces commandes à la main, les unes après les autres, est pénible et source d'erreurs. Et il y a fort à parier que vous aurez complètement oublié ces commandes dans 1 semaine, voire dans 1 heure. Pour autant, c'est parfaitement normal, il n'y a absolument aucun intérêt à se souvenir de toutes ces commandes.

Pour répondre à ces deux problèmes, de gestion de données et d'automatisation, nous allons introduire les notions Bash de variables et de scripts.


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
La commande `echo` affiche une chaîne de caractère, une variable, ou les deux.

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
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script1.sh
```
- Ouvrez le script `script1.sh` avec `nano` et modifiez la variable `sample` avec votre numéro d'échantillon. Sauvegardez le script (`ctrl + o`) et quittez nano (`ctrl + x`).  
Rappel : pas d'espace avant ou après le symbole `=` !
- Lancez le script avec la commande
```
$ bash script1.sh
```

Vérifiez que le déroulement du script se passe bien. Vous avez le temps de prendre un café ☕. Voir plusieurs ☕ 🍪 ☕ 🍪.


## Étape 5 : automatisation de l'analyse : niveau 2

Le script précédent était pratique mais il ne conserve pas les informations liées à l'alignement (nombre de *reads* non-alignés, alignés une fois...).

Proposez une évolution du premier script pour répondre à ce problème. N'hésitez pas à améliorer aussi la lisibilité du script lors de son exécution.

La solution est dans le [script 2](script2.sh). Pour le télécharger, utilisez la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script2.sh
```
Vous remarquerez que la solution proposée pour conserver les informations liées à l'alignement est un peu particulière. Nous allons en discuter, mais dans un premier temps essayer de comprendre l'explication donnée [ici](https://stackoverflow.com/questions/876239/how-can-i-redirect-and-append-both-stdout-and-stderr-to-a-file-with-bash).


## Étape 6 : automatisation de l'analyse : niveau 3 (ninja)

Le script précédent était intéressant mais il ne prend en compte qu'un seul échantillon à la fois. Quel ennui !

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

La ligne `set -euo pipefail` tout au début du script va arrêter celui-ci :
- à la première erreur ;
- si une variable n'est pas définie ;
- si une erreur est rencontrée dans une commande avec un pipe (`|`).

C'est une mesure de sécurité importante pour votre script. Si vous le souhaitez, vous pouvez lire l'article de Aaron Maxwell à ce sujet : [Use the Unofficial Bash Strict Mode (Unless You Looove Debugging)](http://redsymbol.net/articles/unofficial-bash-strict-mode/)

Téléchargez le script 3 avec la commande :
```
$ wget https://raw.githubusercontent.com/omics-school/analyse-rna-seq/master/script3.sh
```

Modifiez le script 3 avec les numéros d'échantillons que vous avez à analyser. Faites bien attention à la variable concernée et sa syntaxe.

Si vous pensez en avoir le temps, lancez le script 3. Comme ce script va automatiser toute l'analyse, il va fonctionner plusieurs minutes et vous aurez peut-être besoin de fermez votre session. Pour ne pas arrêter brutalement l'analyse à la fermture de la session, lancez le script de cette manière :

```
$ nohup bash script3.sh &
```

Le message 
```
nohup: ignoring input and appending output to 'nohup.out'
```
vous rapelle que les messages qui apparaissaient habituellement à l'écran seront redirigés dans le fichier `nohup.out`.


#### Remarque

Au tout début de l'activité, vous avez renommé les fichiers contenant les *reads* (`.fastq.gz`). Cette étape manuelle peut être automatisée, par exemple avec la commande suivante :

```
for name in $(ls 140317*.fastq.gz); do newname=$(echo ${name} | sed 's/140317_SN365_A_L00._//'); mv ${name} ${newname}; done
```

Essayez de comprendre son fonctionnement.
