# ST2 (Théorie des Jeux) - EI Algorithmique Génétique

# Description 

Création de deux programmes pour rendre um plasmide circulaire en utilisant programmation orienté á objet.

# Membres de l'équipe : 
- Oliveira Lima Lawson
- 
- 
- 

# Projet : Codage d'un jeu de sudoku

But du projet : L'objectif est d'utiliser deux approches différentes pour rendre une séquence d'adn circulaire, plus specificment, l'algorithme Monte Carlo Tree Search et algoritmhs genetique.

###########################################################################################################

# Voici les modules qu'il est nécessaire d'importer pour utiliser notre projet

- mathutils
- numpy
- matplotlib

# Utilisation de notre projet

Se placer dans le fichier sudoku_game_by_team6/Sudoku_MVP_lignes_de_commande.py et l'exécuter.

Vous aurez le choix entre reprendre une partie sauvegardée et une nouvelle partie. Vous aurez aussi la possibilité de choisir le niveau de difficulté de la grille, ainsi que les dimensions de celle-ci. Une fois la grille affichée, vous pourrez remplir une case en rentrant sa valeur dans le terminal, ainsi que ses coordonnées. Vous pouvez à tout moment sauvegarder et quitter la partie, ou simplement abandonner. Si vous abandonnez, la solution de la grille s'affichera dans le terminal.


###########################################################################################################
# Organisation des fichiers et modules

- Algo_gen : ce dossier contient les codes du Algorithm Genetique
    - Benchmark : ce fichier containt notre class de benchmark pour l'ago genetique
    - Chromosome : ce fichier containt notre Chromosoome, elle represente un ensemble de 16 gènes
    - Population : ce fichier containt notre Population, elle represente un ensemble de 100 chromosomes
    - Genetics Algo : ce fichier containt notre fonction qu'utilise les classes, elle faire la liason entre les principales fonctions
    - Traj3D : ce fichier containt notre fonction que applique chaque individu de la population (notre chromosome) dans la séquence d'adn
    - Test : ce dossier contient tous les tests qu'on a faits pour Algorithme Genetique

- Solutions_Genetics_Algo : ce dossier contient les solutions obtenus en utilisant un algorithme genetique

- Solutions_Monte_Carlo : ce dossier contient les solutions obtenus en utilisant le MCTS

- Monte_Carlo : ce dossier contient les codes du MCTS

###########################################################################################################

# Nos points fortes
    - Deux programmes très puissance qui obtient des solution optimales
    - Un façon différent de faire le crossover 
    - Des fonctions objectives qui augment la qualité de la solution


###########################################################################################################

# Suggestions pour améliorer notre projet à l'avenir
    - Ajouter une fonction de créer plusieurs populations initialles et ensuite faire un tournoi pour choisir la premier population
    - Ajouter un bouton pour annuler un mouvement
    - Ajouter des modes de jeu pour grilles 16x16 et 25x25

