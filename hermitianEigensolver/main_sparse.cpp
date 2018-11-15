#include <Eigen/Eigenvalues>
#include <Eigen/SparseCore>

typedef Eigen::SparseMatrix<double> Matrix_rds;
typedef Eigen::SparseVector<double> Vector_rds;
typedef Eigen::SparseMatrix<complex<double>> Matrix_cds;
typedef Eigen::SparseVector<complex<double>> Vector_cds;
typedef Eigen::VectorXcd Vector_cdd;

const double DBL_EPSILON = std::numeric_limit<double>::epsilon();

Matrix_rds lanczos_Sparse(Matrix_cds M, int L) {
    int N = M.innerSize();

    // May need to initialise these to zero?
    //Matrix to store vector sequence V
    Matrix_cds V;

    // Tridiagonal system to be diagonalised
    Matrix_rds T;

    // Coefficients in tridiagonal matrix T
    Vector_rds a, b;
    double b0;
    
    // Random initial vector. We still ought to store this densely
    Vector_cdd v = Eigen::VectorXcd::Random(N);
    // Intermediate vectors
    Vecor_rds w,v_gs;
    // Intermediate values
    double temp_norm;

    // Other flags
    is_orthogonal = false;

    b0 = v.norm();
    v /= b0;

    V.col(0)    = v;

    w           = M.selfadjointView<>()*v; // assuming all coefficients of M are stored

    a(0)        = (v.dot(w)).real();
    w           = w - a(0)*v;

    // Lanczos iteration
    for (int i=1; i<L; i++){
        b(i-1) = w.norm();

    }

}
