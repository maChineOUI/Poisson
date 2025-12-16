/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

// void eig_poisson1D(double* eigval, int *la){
//   // TODO: Compute all eigenvalues for the 1D Poisson operator
// }
void eig_poisson1D(double* eigval, int *la){
    int n = *la;
    const double pi = acos(-1.0);
    double theta = pi / (n + 1);

    for (int j = 1; j <= n; ++j) {
        eigval[j - 1] = 2.0 - 2.0 * cos(j * theta);
    }
}


// double eigmax_poisson1D(int *la){
//   // TODO: Compute and return the maximum eigenvalue for the 1D Poisson operator
//   return 0;
// }
double eigmax_poisson1D(int *la){
    int n = *la;
    const double pi = acos(-1.0);
    double theta = pi / (n + 1);
    return 2.0 - 2.0 * cos(n * theta);
}


// double eigmin_poisson1D(int *la){
//   // TODO: Compute and return the minimum eigenvalue for the 1D Poisson operator
//   return 0;
// }
double eigmin_poisson1D(int *la){
    const double pi = acos(-1.0);
    double theta = pi / (*la + 1);
    return 2.0 - 2.0 * cos(1 * theta);
}


// double richardson_alpha_opt(int *la){
//   // TODO: Compute alpha_opt
//   return 0;
// }
double richardson_alpha_opt(int *la){
    double lam_min = eigmin_poisson1D(la);
    double lam_max = eigmax_poisson1D(la);
    return 2.0 / (lam_min + lam_max);
}


/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
// void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
// //   TODO: Implement Richardson iteration
// //   1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
// //   2. Update x = x + alpha*r (use daxpy)
// //   3. Check convergence: ||r||_2 < tol (use dnrm2)
// //   4. Store residual norm in resvec and repeat
// }
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich,
                      int *lab, int *la, int *ku, int*kl,
                      double *tol, int *maxit, double *resvec, int *nbite)
{
    int n = *la;
    int ldA = *lab;
    int inc = 1;
    double one = 1.0, zero = 0.0;
    double alpha = *alpha_rich;

    double *r  = (double*) malloc(n * sizeof(double));
    double *Ax = (double*) malloc(n * sizeof(double));

    double norm_b = dnrm2_(la, RHS, &inc);

    int k;
    for (k = 0; k < *maxit; ++k) {

        /* Ax = A * X */
        dgbmv_("N", la, la, kl, ku, &one, AB, &ldA,
               X, &inc, &zero, Ax, &inc);

        /* r = b - Ax */
        for (int i = 0; i < n; ++i)
            r[i] = RHS[i] - Ax[i];

        double norm_r = dnrm2_(la, r, &inc);
        resvec[k] = norm_r / norm_b;

        if (resvec[k] < *tol) break;

        /* X = X + alpha * r */
        daxpy_(la, &alpha, r, &inc, X, &inc);
    }

    *nbite = k + 1;

    free(r);
    free(Ax);
}


/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
// void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
//   // TODO: Extract diagonal elements from AB and store in MB
//   // MB should contain only the diagonal of A
// }
void extract_MB_jacobi_tridiag(double *AB, double *MB,
                               int *lab, int *la, int *ku, int*kl, int *kv)
{
    int n = *la;
    int ldA = *lab;
    int kv0 = *kv;

    /* 先全部置 0 */
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < ldA; ++i)
            MB[i + j * ldA] = 0.0;

    /* 仅复制主对角线 */
    for (int j = 0; j < n; ++j) {
        int idx = indexABCol(kv0, j, lab);
        MB[idx] = AB[idx];
    }
}


/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
// void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
//   // TODO: Extract diagonal and lower diagonal from AB
//   // MB should contain the lower triangular part (including diagonal) of A
// }
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB,
                                     int *lab, int *la, int *ku, int*kl, int *kv)
{
    int n = *la;
    int ldA = *lab;
    int kv0 = *kv;

    /* 全部置 0 */
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < ldA; ++i)
            MB[i + j * ldA] = 0.0;

    /* 主对角线 */
    for (int j = 0; j < n; ++j) {
        int idx = indexABCol(kv0, j, lab);
        MB[idx] = AB[idx];
    }

    /* 下对角线 (i = j+1) */
    if (*kl >= 1) {
        for (int j = 0; j < n-1; ++j) {
            int idx = indexABCol(kv0 + 1, j, lab);
            MB[idx] = AB[idx];
        }
    }
}


/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
// void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
//   // TODO: Implement Richardson iterative method
// }
void richardson_MB(double *AB, double *RHS, double *X, double *MB,
                   int *lab, int *la, int *ku, int*kl,
                   double *tol, int *maxit, double *resvec, int *nbite)
{
    int n = *la;
    int ldA = *lab;
    int kv = *ku;  /* 主对角线所在行 */
    int inc = 1;
    double one = 1.0, zero = 0.0;

    double *r  = (double*) malloc(n * sizeof(double));
    double *Ax = (double*) malloc(n * sizeof(double));
    double *z  = (double*) malloc(n * sizeof(double));

    double norm_b = dnrm2_(la, RHS, &inc);

    int k;
    for (k = 0; k < *maxit; ++k) {

        /* Ax = A X */
        dgbmv_("N", la, la, kl, ku, &one, AB, &ldA,
               X, &inc, &zero, Ax, &inc);

        /* r = b - Ax */
        for (int i = 0; i < n; ++i)
            r[i] = RHS[i] - Ax[i];

        double norm_r = dnrm2_(la, r, &inc);
        resvec[k] = norm_r / norm_b;

        if (resvec[k] < *tol) break;

        /* Solve z = M^{-1} r */

        /* Jacobi：M 只有对角线 kl=0 */
        if (*kl == 0) {
            for (int i = 0; i < n; ++i) {
                int idx = indexABCol(kv, i, lab);
                z[i] = r[i] / MB[idx];
            }
        }
        else {
            /* GS 前代求解：下三对角 */
            for (int i = 0; i < n; ++i) {
                double sum = r[i];

                if (i > 0) {
                    int idx_sub = indexABCol(kv+1, i-1, lab);
                    sum -= MB[idx_sub] * z[i-1];
                }

                int idx_d = indexABCol(kv, i, lab);
                z[i] = sum / MB[idx_d];
            }
        }

        /* X = X + z */
        daxpy_(la, &one, z, &inc, X, &inc);
    }

    *nbite = k + 1;

    free(r);
    free(Ax);
    free(z);
}


