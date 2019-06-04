class: center, middle

# Analyse de données RNA-seq avec un cluster SLURM

## DU Omiques 2019

Pierre Poulain / @pierrepo

<br /><br /><br /><br /><br /><br />

<div>
<img src="img/logo_DUO.png"
	 height="100px" style="vertical-align:midle;">
 </img>
<div style="display: inline-block; width:100px;"></div>
<img src="img/logo_UPD_USPC.png"
 	 height="120px" style="vertical-align:midle;">
 </img>
</div>

.footer[
Ce contenu est mis à disposition selon les termes de la licence Creative Commons BY-SA 4.0
]

---

layout: true
name: title
class: center, middle
.footer[
DU Omiques 2019
]

---

layout: true
name: contentleft
class: top, left
.footer[
DU Omiques 2019
]

---

layout: true
name: contentcenter
class: top, center
.footer[
DU Omiques 2019
]

---

template: contentleft

# Objectifs de l'activité

--

- Comprendre le fonctionnement du cluster de l'IFB.

- Utiliser les commandes SLURM.

- Adapter un processus d'analyse RNA-seq pour un cluster.

- Lancez un script *sbatch* pour automatiser un processus d'analyse.

- Copier des données entre le serveur et votre ordinateur.

---
template: contentleft

# Dessine-moi un cluster...

--

.center[
<img height="470px" src="img/Penalva__Occigen__Wikimedia__CC-BY-SA.jpg">
]

.footnote[.ref[
	Source : <a href="https://commons.wikimedia.org/wiki/File:Occigen.jpg">Penalva</a>, Wikipedia, CC BY-SA
]]

---
template: contentleft

# Dessine-moi un cluster...

.center[
<img height="520px" src="img/racks_cluster_IFB.jpg">
]

.footnote[.ref[
	Source : Julien Seiler, IFB, CC BY-SA
]]

---
template: contentleft

# Un noeud

.center[
<img height="320px" src="img/node_cluster_IFB.jpg">
]

Un noeud (*node*) = une machine physique = un élément du cluster

Sur ce cluster, un noeud a deux processeurs. Un processeur a 14 coeurs.

.footnote[.ref[
	Source : Julien Seiler, IFB, CC BY-SA
]]


---
template: contentleft

# Some HPC clusters in France

.pure-table.pure-table-bordered.smaller-font[
Cluster | Data center location | Cores | RAM (GB) | Storage (TB) | Access modality
--- | --- | --- | --- | --- | ---
IFB Core | IDRIS - Orsay | 2 000 | 20 008 | 400 | Open to all academic biologists and bioinformaticians
GENOTOUL | Toulouse | 3 064 | 34 304 | 3 000 | Open to all academics with priority to INRA/Occitane region (currently overloaded)
CINES OCCIGEN | Montpellier | 85&nbsp;824 | 202&nbsp;000 | 8 000 | Periodic calls for projects (~2 calls / year)
]

.footnote[.ref[
	Source : Julien Seiler, IFB, CC BY-SA
]]

---
template: contentleft

# SLURM

Simple Linux Utility for Resource Management

---
template: contentleft

# SLURM

.center[
<img height="500px" src="img/slurm_components_IFB.png">
]

.footnote[.ref[
	Source : Julien Seiler, IFB, CC BY-SA
]]

---
template: contentleft

# SLURM

Noeud de soumission de jobs = noeud de connexion

```
$ ssh login@core.cluster.france-bioinformatique.fr
```

--

<br />
<br />
<br />
<br />

On lance un job **depuis** ce noeud mais pas **sur** 😡

s-commands (`srun`, `sbatch`, `squeue`...)

---
template: contentleft

# Doc et liens utiles

[IFB Core Cluster documentation](http://taskforce-nncr.gitlab.cluster.france-bioinformatique.fr/doc/)


---
template: contentleft

# Le processus d'analyse RNA-seq

.center[
<img height="500px" src="img/pipeline_RNA_seq_O_tauri.png">
]

---
template: contentleft

# Conda 🐍

.center[
<img height="400px" src="img/conda.png">
]

---
template: contentleft

# Copie de données 🦄

--

## D
	
<https://filezilla-project.org/>  - logiciel open source

.center[
<img width="1000px" src="img/filezilla.png">
]


Hôte : `sftp://omics-school.net`

Identifiant : `<votre-login-sur-le-serveur>`

Mot de passe : `<votre-mot-de-passe-sur-le-serveur>`

---
template: contentleft

# Copie de données 🦄


## `scp` 

À utiliser dans un *shell* depuis la machine locale.

<br />

De la machine locale vers le serveur : 
```
$ scp fichier.txt ppoulain@omics-school.net:~/repertoire/
```

<br />

Du serveur vers la machine locale : 
```
$ scp ppoulain@omics-school.net:~/repertoire/fichier.txt ./
```


---

template: contentleft

background-color: #cccccc

# C'est parti ! 🚀

## 💻 [Tutoriel](https://omics-school.github.io/analyse-rna-seq/analyse_RNA-seq_O_tauri.html)

## 💻 [Check-list](https://omics-school.github.io/analyse-rna-seq/analyse_RNA-seq_O_tauri_check-list.html)
