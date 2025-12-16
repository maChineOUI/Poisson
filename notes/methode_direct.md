# Problème de Poisson 1D

## Algorithmes mathématiques et implémentation informatique

*(TD/TP de Calcul numérique — Méthodes directes)*

---

## 1. Modèle mathématique (problème continu)

### Énoncé du sujet (français)

> *Soit à résoudre l’équation de la chaleur dans un milieu immobile linéaire homogène avec terme source et isotrope :*
> [
> \begin{cases}
> -k \dfrac{\partial^2 T}{\partial x^2} = g, \quad x \in ]0,1[ \
> T(0) = T_0,; T(1) = T_1
> \end{cases}
> ]

---

### Cas considéré dans le TP

Dans la suite du TP, le cas sans source est retenu :

> *Dans la suite du TP on considère qu’il n’y a pas de source de chaleur, i.e. ( g = 0 ).*

On résout donc :
[
\begin{cases}

* u''(x) = 0,\quad x\in(0,1) \
  u(0)=T_0,\quad u(1)=T_1
  \end{cases}
  ]

La solution analytique est :
[
u(x) = T_0 + x (T_1 - T_0)
]

---

## 2. Discrétisation par différences finies

### Énoncé du sujet

> *On se propose de résoudre cette équation par une méthode de différence finie centrée d’ordre 2.*

---

### Maillage

* Le domaine ([0,1]) est discrétisé en (n+2) points
* Le nombre d’inconnues est :
  [
  n = \texttt{la}
  ]
* Le pas du maillage est :
  [
  h = \frac{1}{n+1}
  ]

Implémentation :

```c
set_grid_points_1D(double* x, int* la);
```

---

### Approximation de la dérivée seconde

Au point intérieur (x_i) :
[
u''(x_i)\approx
\frac{u_{i-1}-2u_i+u_{i+1}}{h^2}
]

L’équation discrète devient :
[
\frac{-u_{i-1}+2u_i-u_{i+1}}{h^2} = 0
]

---

## 3. Système linéaire obtenu

### Énoncé du sujet

> *Pour la suite du TP on notera ce système :*
> [
> A u = f,\quad A\in\mathbb{R}^{n\times n}
> ]

---

### Structure de la matrice

[
A=\frac{1}{h^2}
\begin{pmatrix}
2 & -1 &        \
-1& 2  & -1     \
&\ddots&\ddots&\ddots\
&       & -1 & 2
\end{pmatrix}
]

Propriétés :

* matrice tridiagonale
* symétrique définie positive
* matrice bande

---

### Second membre (conditions de Dirichlet)

Même si (g=0), les conditions aux limites contribuent au second membre :
[
f_1=\frac{T_0}{h^2},\qquad
f_n=\frac{T_1}{h^2}
]

Implémentation :

```c
set_dense_RHS_DBC_1D(double* RHS,
                     int* la,
                     double* BC0,
                     double* BC1);
```

---

## 4. Stockage bande (General Band)

### Énoncé du sujet

> *Une méthode de stockage efficace pour ce type de matrice est le stockage par bande.*

---

### Format General Band (GB)

Pour une matrice tridiagonale :

* (kl = 1) (sous-diagonale)
* (ku = 1) (sur-diagonale)
* (kv = 1) (ligne supplémentaire)
* dimension dominante :
  [
  \texttt{lab} = kv + kl + ku + 1 = 4
  ]

La matrice est stockée dans un tableau :

```c
AB[lab][n]   // column-major
```

Accès aux coefficients :

```c
int indexABCol(int i, int j, int *lab);
```

---

### Construction de l’opérateur de Poisson

```c
set_GB_operator_colMajor_poisson1D(double* AB,
                                   int *lab,
                                   int *la,
                                   int *kv);
```

Cette fonction construit exactement :
[
A=\frac{1}{h^2}\operatorname{tridiag}(-1,2,-1)
]

---

## 5. Méthodes directes de résolution

Les exercices 4 à 6 consistent à **implémenter et comparer trois méthodes directes** pour résoudre le système linéaire.

---

### Méthode 1 : LU bande générique (TRF)

Implémentation LAPACK :

```c
dgbtrf_(&n, &n, &kl, &ku, AB, &lab, ipiv, &info);
dgbtrs_("N", &n, &kl, &ku, &NRHS,
         AB, &lab, ipiv, RHS, &n, &info);
```

Caractéristiques :

* méthode générale
* complexité en temps : (O(n))
* constante de coût relativement élevée

---

### Méthode 2 : LU tridiagonale spécialisée (TRI)

Implémentation spécifique :

```c
dgbtrftridiag(int *la, int *n,
              int *kl, int *ku,
              double *AB, int *lab,
              int *ipiv, int *info);
```

Principe :
[
\begin{aligned}
L_{i+1,i} &= \frac{A_{i+1,i}}{A_{i,i}} \
A_{i+1,i+1} &\leftarrow A_{i+1,i+1}
- L_{i+1,i} A_{i,i+1}
\end{aligned}
]

Caractéristiques :

* exploitant strictement la structure tridiagonale
* complexité (O(n))
* constantes minimales (méthode la plus rapide)

---

### Méthode 3 : Solveur LAPACK tout-en-un (SV)

Implémentation :

```c
dgbsv_(&n, &kl, &ku, &NRHS,
       AB, &lab, ipiv,
       RHS, &n, &info);
```

Remarque :
[
\texttt{dgbsv} = \texttt{dgbtrf} + \texttt{dgbtrs}
]

Caractéristiques :

* interface la plus simple
* performances proches de la méthode TRF

---

## 6. Évaluation de l’erreur

L’erreur est mesurée par l’erreur relative avant :
[
\frac{|u_h-u_{\text{exact}}|*2}
{|u*{\text{exact}}|_2}
]

Implémentation :

```c
double relative_forward_error(double* x,
                              double* y,
                              int* la);
```

---

## 7. Comparaison des méthodes

| Méthode | Implémentation    | Complexité | Remarque         |
| ------- | ----------------- | ---------- | ---------------- |
| TRF     | `dgbtrf + dgbtrs` | (O(n))     | Générale         |
| TRI     | `dgbtrftridiag`   | (O(n))     | Optimale         |
| SV      | `dgbsv`           | (O(n))     | Interface simple |

---




# Poisson 1D 问题：数学算法与计算机实现笔记

*(TD/TP Calcul numérique — Méthodes directes)*

---

## 1. 数学问题描述（连续模型）

### 题目原文

> *Soit à résoudre l’équation de la chaleur dans un milieu immobile linéaire homogène avec terme source et isotrope :*
> [
> \begin{cases}
> -k \dfrac{\partial^2 T}{\partial x^2} = g, \quad x \in ]0,1[ \
> T(0) = T_0,; T(1) = T_1
> \end{cases}
> ]

### 中文对照

在区间 ( (0,1) ) 上求解一维稳态热传导（Poisson）方程，给定 Dirichlet 边界条件。

---

### 本 TP 采用的**特例（与最终代码一致）**

题目后续明确采用：

> *Dans la suite du TP on considère qu’il n’y a pas de source de chaleur, i.e. (g=0).*

即：
[
\begin{cases}

* u''(x) = 0,\quad x\in(0,1) \
  u(0)=T_0,\quad u(1)=T_1
  \end{cases}
  ]

解析解为：
[
u(x)=T_0 + x(T_1-T_0)
]

---

## 2. 数值离散（有限差分）

### 题目原文

> *On se propose de résoudre cette équation par une méthode de différence finie centrée d’ordre 2.*

### 中文对照

采用**二阶中心有限差分格式**。

---

### 网格离散

* 将区间 ([0,1]) 均匀划分为 (n+2) 个点
* 内部未知点数：
  [
  n = \text{la}
  ]
* 网格步长：
  [
  h=\frac{1}{n+1}
  ]

代码实现函数：

```c
set_grid_points_1D(double* x, int* la);
```

---

### 二阶中心差分

在内部点 (x_i)：
[
u''(x_i)\approx\frac{u_{i-1}-2u_i+u_{i+1}}{h^2}
]

代入方程 (-u''=0)，得到离散线性系统：
[
\frac{-u_{i-1}+2u_i-u_{i+1}}{h^2}=0
]

---

## 3. 线性代数形式

### 题目原文

> *Pour la suite du TP on notera ce système :*
> [
> A u = f,\quad A\in\mathbb{R}^{n\times n}
> ]

---

### 矩阵结构

[
A=\frac{1}{h^2}
\begin{pmatrix}
2 & -1 &   &   \
-1& 2  & -1&   \
&\ddots&\ddots&\ddots\
&   & -1 & 2
\end{pmatrix}
]

性质：

* 三对角矩阵
* 对称正定（SPD）
* 带状结构（band matrix）

---

### 右端项（Dirichlet 边界贡献）

尽管 (g=0)，边界条件仍贡献 RHS：

[
f_1=\frac{T_0}{h^2},\qquad
f_n=\frac{T_1}{h^2}
]

代码实现：

```c
set_dense_RHS_DBC_1D(double* RHS, int* la,
                     double* BC0, double* BC1);
```

---

## 4. 带状矩阵存储（GB Storage）

### 题目原文（法语）

> *Une méthode de stockage efficace pour ce type de matrice est le stockage par bande.*

---

### General Band（GB）格式

对于三对角矩阵：

* (kl=1)：下对角线数
* (ku=1)：上对角线数
* (kv=1)：额外存储行
* leading dimension：
  [
  lab = kv + kl + ku + 1 = 4
  ]

矩阵存储为：

```c
AB[lab][n]   // column-major
```

索引函数：

```c
int indexABCol(int i, int j, int *lab);
```

---

### Poisson 1D 矩阵构造

```c
set_GB_operator_colMajor_poisson1D(double* AB,
                                   int *lab,
                                   int *la,
                                   int *kv);
```

实现的正是：
[
A=\frac{1}{h^2}\operatorname{tridiag}(-1,2,-1)
]

---

## 5. 三种直接求解方法（核心）

> **Exercices 4–6 的目标：构建并比较三种直接法的性能**

---

### 方法一：通用带状 LU（TRF）

使用 LAPACK：

```c
dgbtrf_(&n, &n, &kl, &ku, AB, &lab, ipiv, &info);
dgbtrs_("N", &n, &kl, &ku, &NRHS,
         AB, &lab, ipiv, RHS, &n, &info);
```

特点：

* 通用性强
* 时间复杂度 (O(n))
* 常数较大

---

### 方法二：自实现三对角 LU（TRI）

代码函数：

```c
dgbtrftridiag(int *la, int *n, int *kl, int *ku,
              double *AB, int *lab,
              int *ipiv, int *info);
```

算法思想（Thomas-like）：
[
\begin{aligned}
L_{i+1,i} &= \frac{A_{i+1,i}}{A_{i,i}} \
A_{i+1,i+1} &\leftarrow A_{i+1,i+1}
- L_{i+1,i}A_{i,i+1}
\end{aligned}
]

特点：

* 专用于三对角
* 时间复杂度 (O(n))
* 常数最小（最快）

---

### 方法三：一体化 LAPACK（SV）

使用：

```c
dgbsv_(&n, &kl, &ku, &NRHS,
       AB, &lab, ipiv,
       RHS, &n, &info);
```

等价于：
[
\text{dgbsv} = \text{dgbtrf} + \text{dgbtrs}
]

特点：

* 接口最简单
* 性能与 TRF 接近

---

## 6. 误差评估

数值解与解析解比较：

[
\text{relative forward error}
= \frac{|u_h-u_{\text{exact}}|*2}{|u*{\text{exact}}|_2}
]

代码函数：

```c
double relative_forward_error(double* x,
                              double* y,
                              int* la);
```

---

## 7. 总结

| 方法  | 实现函数              | 复杂度    | 特点    |
| --- | ----------------- | ------ | ----- |
| TRF | `dgbtrf + dgbtrs` | (O(n)) | 通用、稳健 |
| TRI | `dgbtrftridiag`   | (O(n)) | 最快、特化 |
| SV  | `dgbsv`           | (O(n)) | 最简接口  |

---

