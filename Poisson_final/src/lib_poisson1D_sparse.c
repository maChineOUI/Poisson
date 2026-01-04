// lib_poisson1D_sparse.c                     
// Sparse storage for 1D Poisson operator
// Format creux CSR / CSC pour Poisson 1D
// 一维 Poisson 算子的稀疏存储（CSR / CSC)

#include <stdlib.h>
#include "lib_poisson1D_sparse.h"


/* 
  CSR storage for 1D Poisson matrix 
  Stockage CSR de la matrice de Poisson 1D 
  一维 Poisson 矩阵的 CSR（Compressed Sparse Row）存储      
*/

/*
  A =
   2 -1  0 ...
  -1  2 -1 ...
   0 -1  2 ...
 
  nnz = 3*n - 2
  Nombre de coefficients non nuls
  非零元素个数
*/

void poisson1D_CSR(int n, double **val, int **col, int **rowptr)
{
  int nnz = 3*n - 2;   // Nombre total de valeurs non nulles
                       // 非零元素总数

  *val    = (double*) malloc(nnz * sizeof(double)); // valeurs
  *col    = (int*)    malloc(nnz * sizeof(int));    // indices de colonnes
  *rowptr = (int*)    malloc((n+1) * sizeof(int));  // pointeurs de lignes
                                                    // 行指针 

  int k = 0;                 // index global dans val/col
                             // val/col 的全局索引
  (*rowptr)[0] = 0;          // début de la première ligne
                             // 第一行起始位置

  for (int i = 0; i < n; i++) {

    // subdiagonal : A(i,i-1) = -1 
    // sous-diagonale 
    // 次对角线 
    if (i > 0) {
      (*val)[k] = -1.0;
      (*col)[k] = i - 1;
      k++;
    }

    // diagonal : A(i,i) = 2 
    // diagonale principale 
    // 主对角线 
    (*val)[k] = 2.0;
    (*col)[k] = i;
    k++;

    // superdiagonal : A(i,i+1) = -1 
    // sur-diagonale 
    // 上对角线 
    if (i < n-1) {
      (*val)[k] = -1.0;
      (*col)[k] = i + 1;
      k++;
    }

    // fin de la ligne i 
    // ligne i terminée 
    // 第 i 行结束 
    (*rowptr)[i+1] = k;
  }
}

// CSR matrix-vector product                                 
// Produit matrice-vecteur en format CSR                     
// CSR 格式下的矩阵-向量乘法                                 
// y = A * x                                                 

void dcsrmv(int n, double *val, int *col, int *rowptr,
            double *x, double *y)
{
  for (int i = 0; i < n; i++) {
    y[i] = 0.0;  // initialisation de y 
                 // 初始化 y 
    for (int k = rowptr[i]; k < rowptr[i+1]; k++) {
      y[i] += val[k] * x[col[k]];
      // accumulation y_i += A_ij * x_j 
      // 累加 y_i += A_ij * x_j 
    }
  }
}

// CSC storage for 1D Poisson matrix                         
// Stockage CSC de la matrice de Poisson 1D                  
// 一维 Poisson 矩阵的 CSC（Compressed Sparse Column）存储   

void poisson1D_CSC(int n, double **val, int **row, int **colptr)
{
  int nnz = 3*n - 2;   // nombre de coefficients non nuls 
                       // 非零元素个数 

  *val    = (double*) malloc(nnz * sizeof(double)); // valeurs 
  *row    = (int*)    malloc(nnz * sizeof(int));    // indices de lignes 
  *colptr = (int*)    malloc((n+1) * sizeof(int));  // pointeurs de colonnes 
                                                    // 列指针 

  int k = 0;                 // index global 
                             // 全局索引 
  (*colptr)[0] = 0;          // début de la première colonne 
                             // 第一列起始 

  for (int j = 0; j < n; j++) {

    // subdiagonal : A(j-1,j) = -1 
    // sous-diagonale 
    // 次对角线 
    if (j > 0) {
      (*val)[k] = -1.0;
      (*row)[k] = j - 1;
      k++;
    }

    // diagonal : A(j,j) = 2 
    // diagonale principale 
    // 主对角线 
    (*val)[k] = 2.0;
    (*row)[k] = j;
    k++;

    // superdiagonal : A(j+1,j) = -1 
    // sur-diagonale 
    // 上对角线 
    if (j < n-1) {
      (*val)[k] = -1.0;
      (*row)[k] = j + 1;
      k++;
    }

    // fin de la colonne j 
    // colonne j terminée 
    // 第 j 列结束 
    (*colptr)[j+1] = k;
  }
}

// CSC matrix-vector product                                 
// Produit matrice-vecteur en format CSC                     
// CSC 格式下的矩阵-向量乘法                                 
// y = A * x                                                 

void dcscmv(int n, double *val, int *row, int *colptr,
            double *x, double *y)
{
  // initialisation du vecteur résultat 
  // 初始化结果向量 
  for (int i = 0; i < n; i++)
    y[i] = 0.0;

  for (int j = 0; j < n; j++) {
    for (int k = colptr[j]; k < colptr[j+1]; k++) {
      y[row[k]] += val[k] * x[j];
      // accumulation y_i += A_ij * x_j 
      // 累加 y_i += A_ij * x_j 
    }
  }
}
