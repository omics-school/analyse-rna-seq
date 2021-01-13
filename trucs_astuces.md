# Trucs et astuces

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

