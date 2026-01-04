/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

// void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
//   // TODO: Fill AB with the tridiagonal Poisson operator
// }

// Build Poisson 1D operator in GB column-major format
// 构造一维 Poisson 方程的 GB（带状）列主序矩阵
// A = (1/h^2) * tridiag(-1, 2, -1)

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv)
{
    int n    = *la;        // Number of unknowns / 未知数个数
    int ldab = *lab;       // Leading dimension of AB / AB 的 leading dimension

    double h    = 1.0 / (n + 1);        // Grid spacing / 网格步长
    double diag =  2.0 / (h * h);       // Main diagonal value / 主对角线
    double off  = -1.0 / (h * h);       // Off-diagonal value / 副对角线

    // Initialize AB to zero
    // 初始化带状矩阵为 0
    for (int j = 0; j < n; j++)
        for (int i = 0; i < ldab; i++)
            AB[indexABCol(i, j, lab)] = 0.0;

    // Explicit GB row indices for tridiagonal matrix
    // 三对角矩阵在 GB 存储中的行索引
    int row_super = 1;   // upper diagonal / 上对角
    int row_diag  = 2;   // main diagonal  / 主对角
    int row_sub   = 3;   // lower diagonal / 下对角

    // Fill the tridiagonal Poisson operator
    // 填充 Poisson 三对角算子
    for (int j = 0; j < n; j++) {

        if (j > 0)
            AB[indexABCol(row_super, j, lab)] = off;

        AB[indexABCol(row_diag, j, lab)] = diag;

        if (j < n - 1)
            AB[indexABCol(row_sub, j, lab)] = off;
    }
}


// void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
//   // TODO: Fill AB with the identity matrix
//   // Only the main diagonal should have 1, all other entries are 0
// }

// Build identity matrix in GB format
// 构造 GB 格式下的单位矩阵

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv)
{
    for (int j = 0; j < *la; j++) {
        for (int i = 0; i < *lab; i++)
            AB[indexABCol(i, j, lab)] = 0.0;

        // Set main diagonal to 1
        // 主对角线设为 1
        AB[indexABCol(*kv, j, lab)] = 1.0;
    }
}


// void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
//   // TODO: Compute RHS vector
// }  

// Build RHS for Poisson 1D with Dirichlet BC
// 构造一维 Poisson 方程（Dirichlet 边界）的右端项
// Model: -u'' = 0

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1)
{
    double h  = 1.0 / ((*la) + 1);   // Grid spacing / 网格步长
    double h2 = h * h;

    // Initialize RHS to zero
    // RHS 初始化为 0
    for (int i = 0; i < *la; i++)
        RHS[i] = 0.0;

    // Boundary contributions
    // 边界条件对 RHS 的贡献
    RHS[0]     += (*BC0) / h2;
    RHS[*la-1] += (*BC1) / h2;
}


// void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
//   // TODO: Compute the exact analytical solution at each grid point
//   // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
// }  

// Analytical solution for Poisson 1D with g = 0
// 一维 Poisson 方程（无源项）的解析解
// u(x) = BC0 + x (BC1 - BC0)

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X,
                                    int* la, double* BC0, double* BC1)
{
    for (int i = 0; i < *la; i++)
        EX_SOL[i] = (*BC0) + X[i] * ((*BC1) - (*BC0));
}


// void set_grid_points_1D(double* x, int* la){
//   // TODO: Generate uniformly spaced grid points in [0,1]
// }

// Generate uniformly spaced grid points in (0,1)
// 在区间 (0,1) 上生成均匀网格点

void set_grid_points_1D(double* x, int* la)
{
    double h = 1.0 / ((*la) + 1);
    for (int i = 0; i < *la; i++)
        x[i] = (i + 1) * h;
}

// double relative_forward_error(double* x, double* y, int* la){
//   // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
//   return 0.0;
// }

// Relative forward error: ||x - y|| / ||y||
// 相对前向误差：||x - y|| / ||y||

double relative_forward_error(double* x, double* y, int* la)
{
    double num = 0.0, den = 0.0;

    for (int i = 0; i < *la; i++) {
        num += (x[i] - y[i]) * (x[i] - y[i]);
        den += y[i] * y[i];
    }

    return sqrt(num) / sqrt(den);
}

// int indexABCol(int i, int j, int *lab){
//   // TODO: Return the correct index formula for column-major band storage
//   return 0;
// }

// Column-major GB index
// GB 列主序存储的索引计算

int indexABCol(int i, int j, int *lab)
{
    return j * (*lab) + i;
}

// int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
//   // TODO: Implement specialized LU factorization for tridiagonal matrices
//   return *info;
// }

// Specialized LU factorization for tridiagonal matrices
// 三对角矩阵的专用 LU 分解（GB 格式，O(n)）
// No pivoting (Poisson 1D is SPD)
// 不使用主元选取（Poisson 1D 为对称正定）

int dgbtrftridiag(int *la, int *n, int *kl, int *ku,
                  double *AB, int *lab, int *ipiv, int *info){
  int N    = *n;
  int ldab = *lab;
  int KL   = *kl;
  int KU   = *ku;
  (void)ldab;  

  if (KL != 1 || KU != 1) {
    *info = -1;
    return *info;
  }


  int row_super = 1;  
  int row_diag  = 2;  // A(i,   i) / U(i,   i)
  int row_sub   = 3;  // A(i+1, i) / L(i+1, i)

  for (int i = 0; i < N; i++){
    ipiv[i] = i + 1;
  }

    /* 
    * LU factorization for tridiagonal matrix in LAPACK GB format
    * 三对角矩阵的 LU 分解（LAPACK 带状存储，原地进行）
    */

    /* Assume success by default (LAPACK convention)
    * 默认假设分解成功（LAPACK 约定）
    */
    *info = 0;

    /* Loop over columns j = 0 ... N-2
    * 对每一列做一次消元（最后一行无需消元）
    */
    for (int j = 0; j < N - 1; j++) {

        /* Pointer to pivot U(j,j)
        * 指向当前主元 U(j,j)
        */
        double *diag_j_ptr = &AB[indexABCol(row_diag, j, lab)];
        double  diag_j     = *diag_j_ptr;

        /* Check for zero (or nearly zero) pivot
        * 主元检查，避免除零
        */
        if (fabs(diag_j) < DBL_EPSILON) {
            *info = j + 1;   // 1-based index (LAPACK)
            return *info;
        }

        /* Compute elimination factor:
        * L(j+1,j) = A(j+1,j) / U(j,j)
        * 计算消元因子，并原地存入下对角
        */
        double *sub_j_ptr = &AB[indexABCol(row_sub, j, lab)];
        double  l         = (*sub_j_ptr) / diag_j;
        *sub_j_ptr        = l;

        /* Extract superdiagonal U(j,j+1)
        * 读取上对角元素 U(j,j+1)
        * (stored in column j+1 in GB format)
        */
        double super_j = AB[indexABCol(row_super, j + 1, lab)];

        /* Update next pivot:
        * U(j+1,j+1) = A(j+1,j+1) - L(j+1,j)*U(j,j+1)
        * 更新下一主对角（三对角 LU 的核心递推）
        */
        double *diag_j1_ptr = &AB[indexABCol(row_diag, j + 1, lab)];
        *diag_j1_ptr -= l * super_j;
    }

    /* Check last pivot U(N-1,N-1)
    * 检查最后一个主元
    */
    double last_diag = AB[indexABCol(row_diag, N - 1, lab)];
    if (fabs(last_diag) < DBL_EPSILON) {
        *info = N;
    }

    /* Return LAPACK-style status code
    * 返回 LAPACK 风格状态码
    */
    return *info;
}