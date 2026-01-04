# Projet Poisson 1D — Méthodes directes et itératives

Ce dépôt contient un projet de résolution du **problème de Poisson 1D** (différences finies), avec comparaison entre une **méthode directe** (LU via LAPACK) et des **méthodes itératives** (Richardson, Jacobi, Gauss–Seidel).

## Organisation du dépôt

- **`Poisson_seance2/`**  
  Contient le travail correspondant aux **six premiers exercices** (séance 2), servant de base progressive (mise en place, premières fonctions, tests/notes, etc.).

- **`Poisson_final/`**  
  Contient la **version finale du projet**.

  - **`Poisson_final/include/`** et **`Poisson_final/src/`** : **code complet** du projet (bibliothèque + exécutables).
  - **`Poisson_final/results/`** : résultats obtenus lors de l’exécution sur ma machine, incluant notamment les **deux figures de convergence** (`convergence_resvec.png` et `convergence_resvec.pdf`).
  - **`Poisson_final/rapport.md`** : **rapport du projet** (description, méthodes, protocole, résultats et analyse).

## Exécution (version finale)
Les scripts et fichiers nécessaires à l’exécution se trouvent dans `Poisson_final/` (Makefile, Dockerfile, scripts).  
Voir `Poisson_final/README` pour les commandes de compilation/exécution.
