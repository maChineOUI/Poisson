## Rapport de projet — Résolution du problème de Poisson 1D

**Méthodes directes et itératives : LU (LAPACK), Richardson, Richardson préconditionné (Jacobi, Gauss–Seidel)**

### 1. Objet du projet

Ce projet vise à résoudre un problème de Poisson 1D (modèle stationnaire de diffusion / chaleur) discrétisé par différences finies. L’étude compare :

* une **méthode directe** via factorisation **LU** (LAPACK) ;
* des **méthodes itératives** de type **Richardson** :

  1. Richardson avec paramètre (\alpha) choisi par une formule spectrale ((\alpha_{\text{opt}})),
  2. Richardson préconditionné **Jacobi**,
  3. Richardson préconditionné **Gauss–Seidel**.

L’implémentation impose l’usage d’un stockage **bande (General Band / GB)**, la production de sorties numériques (solution, résidus) et une visualisation de la convergence.

---

### 2. Modèle et discrétisation

On considère :
[
-u''(x)=f(x)\quad \text{sur } [0,1],\qquad u(0)=T_0,; u(1)=T_1.
]
En discrétisant sur une grille uniforme et en ne conservant que les inconnues **intérieures**, on obtient un système linéaire :
[
A x = b,
]
où (A) est **tridiagonale** (cas Poisson 1D standard) :

* diagonale : (2),
* sous- et sur-diagonales : (-1).

Dans le cas (f(x)=0), la solution continue compatible avec les seules conditions de bord est linéaire :
[
u(x)=T_0 + x,(T_1-T_0),
]
et elle est évaluée dans le code (fichier `EX_SOL`) comme référence qualitative.

---

### 3. Organisation du projet et reproductibilité

**Arborescence :**

* `src/` : solveurs et bibliothèque C
* `include/` : en-têtes (BLAS/LAPACK via headers ATLAS)
* `docker/Dockerfile` : environnement reproductible (Ubuntu 20.04)
* `results/` : sorties (`AB.dat`, `RHS.dat`, `SOL.dat`, `RESVEC.dat`, etc.)
* `bin/` : exécutables (`tpPoisson1D_direct`, `tpPoisson1D_iter`)

**Environnement Docker :**

* Ubuntu 20.04, `gcc 9.4.0`
* bibliothèques BLAS/LAPACK dynamiques (ex. `libblas`, `liblapack`)

---

### 4. Implémentation numérique

#### 4.1 Stockage GB (General Band)

La matrice est stockée en **col-major** au format bande, avec :
[
ku=1,\quad kl=1,\quad kv=0,\quad lab = kv+kl+ku+1=3.
]
La multiplication matrice-vecteur (Ax) est réalisée par :

* `cblas_dgbmv`.

---

#### 4.2 Méthode de Richardson (pas fixe)

Schéma :
[
x^{(k+1)} = x^{(k)} + \alpha,(b-Ax^{(k)}),\qquad r^{(k)} = b-Ax^{(k)}.
]
Critère d’arrêt (résidu relatif) :
[
\frac{|r^{(k)}|_2}{|b|_2} < \texttt{tol}.
]

Routines BLAS utilisées :

* `dgbmv` : calcul de (Ax),
* `dnrm2` : normes (|r|_2) et (|b|_2),
* `daxpy` : mise à jour (x \leftarrow x + \alpha r).

---

#### 4.3 Choix de (\alpha) par estimation spectrale

Pour l’opérateur de Poisson 1D, les valeurs propres discrètes (points intérieurs) satisfont :
[
\lambda_k = 2 - 2\cos!\left(\frac{k\pi}{n+1}\right),\quad k=1,\dots,n.
]
On en déduit :
[
\alpha_{\text{opt}}=\frac{2}{\lambda_{\min}+\lambda_{\max}}.
]
Dans notre configuration (grand (n)), on a typiquement (\lambda_{\max}\approx 4) et (\lambda_{\min}) très petit (de l’ordre de (\pi^2/(n+1)^2)), ce qui conduit à :
[
\alpha_{\text{opt}} \approx \frac{2}{0+4} = 0.5
]
(approximation asymptotique, suffisante pour interpréter les résultats). Le programme affiche d’ailleurs (\alpha_{\text{opt}}\approx 0.5).

---

#### 4.4 Richardson préconditionné (Jacobi / Gauss–Seidel)

Forme générale :
[
x^{(k+1)} = x^{(k)} + M^{-1}(b-Ax^{(k)}).
]

* **Jacobi** : (M=D), la diagonale de (A).
* **Gauss–Seidel** : (M=D-E) (diagonale + partie strictement inférieure).

L’extraction de (M) en GB est faite par :

* `extract_MB_jacobi_tridiag`,
* `extract_MB_gauss_seidel_tridiag`.

Application :

* Jacobi : (z_i = r_i / D_{ii}),
* GS : résolution triangulaire par substitution avant ((D-E)z=r), puis (x \leftarrow x + z).

---

### 5. Ajustements effectués (fiabilité et traçabilité)

Plusieurs corrections/choix ont été nécessaires pour obtenir des résultats cohérents :

1. **Paramètres unifiés** entre direct et itératif (même problème numérique)
   `nbpoints = 10000`, `T0 = 5`, `T1 = 20`.

2. **Cohérence du stockage bande** pour le préconditionneur `MB`
   Correction : `kv=0`, `ku=1`, `kl=1`, `lab=3` également pour `MB`.

3. **Robustesse compilation**
   Ajout de `#include <string.h>` (usage de `memset`).

4. **Reporting expérimental**
   Affichage : `Method`, `nbite`, `final relres`, `tol`, `maxit`, et sauvegarde de l’historique `RESVEC.dat` pour tracer la convergence.

---

### 6. Protocole expérimental

Paramètres principaux :

* `nbpoints = 10000` (\Rightarrow n = la = nbpoints - 2),
* `tol = 1e-3`, `maxit = 5000`,
* initialisation : `SOL = 0`.

Cas testés côté itératif :

* (0) Richardson avec (\alpha_{\text{opt}}),
* (1) Richardson préconditionné Jacobi,
* (2) Richardson préconditionné Gauss–Seidel.

Côté direct :

* solveur LU LAPACK.

Sorties : `SOL*.dat`, `EX_SOL*.dat`, `RHS*.dat`, `RESVEC_*.dat` + figure `convergence_resvec.png`.

---

### 7. Résultats

#### 7.1 Méthode directe (référence)

Exécution `tpPoisson1D_direct 0` (LU LAPACK) : erreur/indicateur relatif reporté très faible (ordre (10^{-11})–(10^{-12})), ce qui valide l’assemblage (matrice GB, RHS et conditions de Dirichlet).

#### 7.2 Méthodes itératives (résidu relatif)

Pour `tol = 1e-3`, `maxit = 5000` :

* **Richardson (\alpha_{\text{opt}})** :
  `nbite = 5000`, `final relres ≈ 1.263285e-03` (ne passe pas sous la tolérance).

* **Jacobi** :
  `nbite = 5000`, `final relres ≈ 1.263285e-03` (même comportement).

* **Gauss–Seidel** :
  `nbite = 2174`, `final relres ≈ 9.998630e-04` (atteint le critère).

---

### 8. Pourquoi (\alpha_{\text{opt}}) et Jacobi donnent la même courbe ?

Ici, la matrice de Poisson 1D a une diagonale constante :
[
D = 2I.
]
Le préconditionneur Jacobi étant (M=D), on a :
[
M^{-1} = \frac{1}{2}I.
]
Donc Richardson préconditionné Jacobi devient :
[
x^{(k+1)} = x^{(k)} + \frac{1}{2}(b-Ax^{(k)}),
]
c’est-à-dire **Richardson avec (\alpha=0.5)**.

Or, dans notre cas (grand (n)), l’estimation spectrale donne aussi (\alpha_{\text{opt}} \approx 0.5). Les deux méthodes réalisent donc **la même mise à jour**, d’où la superposition parfaite des courbes de résidu.

---

### 9. Pourquoi Gauss–Seidel converge plus vite ?

Gauss–Seidel intègre une partie de la structure de (A) via un préconditionneur triangulaire ((D-E)). Pour le Poisson 1D, cette prise en compte réduit le facteur de convergence par rapport à Jacobi.

Un point important est que, lorsque (n) est grand, (\lambda_{\min}) devient très petit : le conditionnement (au sens spectral) se dégrade fortement, ce qui rend Richardson/Jacobi **très lent** (même si (\alpha) est “optimal” au sens classique). Gauss–Seidel améliore la situation et atteint ici `tol` en 2174 itérations.

---

### 10. Éléments à inclure dans le rapport

* Figure principale : `results/convergence_resvec.png` (échelle log en Y + ligne `tol`).
* Optionnel : comparaison `SOL.dat` vs `EX_SOL.dat` (ou vs solution directe) pour illustrer la qualité de solution.

---

### 11. Conclusion

Le projet met en évidence deux points :

1. **Correction** : le solveur direct LU fournit une référence numérique très précise, validant l’assemblage et le format bande.
2. **Comportements attendus** :

* Jacobi et Richardson((\alpha=0.5)) coïncident pour cette matrice (diagonale (2I)), ce qui explique la superposition des courbes ;
* Gauss–Seidel converge plus rapidement et atteint la tolérance dans le budget d’itérations.

Ces résultats sont cohérents avec l’analyse spectrale et avec les propriétés classiques des schémas itératifs sur le Poisson 1D.

---