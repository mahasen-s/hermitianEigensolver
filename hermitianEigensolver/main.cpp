#include <Eigen/Eigenvalues>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <algorithm> // needed for min,
#include <random>
#include <iostream>
#include <iomanip>
#include <complex>
#include <cmath>
#include <chrono>
#include <tuple>
//#include <type_traits>

//using namespace std::complex_literals;

typedef Eigen::MatrixXcd Matrix_cdd; 	// cdd = complex, double, dense
typedef Eigen::VectorXcd Vector_cdd;
typedef Eigen::MatrixXd Matrix_rdd; 	// rdd = real, double, dense
typedef Eigen::VectorXd Vector_rdd;
typedef Eigen::Matrix<bool,Eigen::Dynamic,1> Vector_bool;

typedef Eigen::SparseMatrix<std::complex<double>> 	Matrix_cds; 	// cds = complex, double, sparse
typedef Eigen::SparseVector<std::complex<double>> 	Vector_cds;
typedef Eigen::SparseMatrix<double> Matrix_rds; 	// cds = complex, double, sparse
typedef Eigen::SparseVector<double> Vector_rds;


const double DBL_EPSILON =std::numeric_limits<double>::epsilon();

template <class C, class D>
C projector(const D &V, const C & r, long ind) {
    // Project r into V
    // TODO: test that this works!
    long N  = V.innerSize();
    C rOut  = (V.block(0,0,N,ind)*((V.block(0,0,N,ind)).transpose()))*r;
    return rOut;
}

template <class C, class D>
void projector2(const D &V, const C & r, long ind, C & rOut, C & bOut) {
    // Project r into V and compute V'*r
    // TODO: test that this works!
    long N  = V.innerSize();
    rOut  = projector(V, r, ind);
    bOut  = ((V.block(0,0,N,ind)).transpose())*r;
    return;
}

//template <class C, class D>
void robustOrthogonalise(const Matrix_cdd &V, Vector_cdd &r,long index, double &normRes, bool &stopAlg) {
    // C and D should both be both real or both complex floating precision.
    // C should be a matrix, D should be a vector
    Vector_cdd r_temp;
    Vector_cdd b_temp;
    double normr0   = r.norm();
    double rMr;
    int numReorths  = 0;
    long N          = r.innerSize();

    // Reorthogonalise
    r               = r - projector(V,r,index);
    normRes         = r.norm();
    numReorths++;

    while( (normRes<=(normr0/sqrt(2.0)) && numReorths<5)){
        r           = r - projector(V,r,index);
        normr0      = normRes;
        normRes     = r.norm();
        numReorths++;
    }

    if (normRes<=(normr0/sqrt(2.0))) {
        // Cannot reorthogonalise, invariant subspace found
        normRes     = 0;
        stopAlg     = true;

        for (int restart = 0; restart < 3; ++restart) {
            // Do a random restart, try at most 3 times
            r       = Eigen::VectorXcd::Random(N);

            // Orthongalise r
            r       = r - projector(V,r,index);
            rMr     = r.norm();
            r       /=rMr;

            // Reorthonalise if necessary
            for (int reorth = 0; reorth < 5; ++reorth) {
                // Check orthogonalisty
                projector2(V,r,index,r_temp,b_temp);
                rMr                         = r_temp.norm();

                if ( (abs(rMr-1) <= 1.0e-10) && ((b_temp.cwiseAbs().array() <= 1.0e-10).all()) ) {
                    stopAlg     = true;
                    return;
                }

                // Re orthogonalise
                r           = r - r_temp;
                r           /= r.norm();

            }
        }
    } else {
        r   /= normRes;
    }
};

struct krylovSchurReturnObject {
    Vector_rdd eigenvalues;
    Matrix_rdd eigenvectors;
    Vector_bool isConverged;
};

template <class C>
krylovSchurReturnObject krylovSchur(const C &M, const long k0, const long maxIter = 300, const double tol = 1.0e-14) {
    // M is the Hermitian matrix we want to diagonalise
    // M must inherit from Matrix_cdd or Matrix_cds
    // k is the number of eigenvalues/eigenvectors we want to find
    // maxIter is the maximum number of iterations allowed
    // tol is the numerical tolerance of solutions allowed
    // TODO: Started using longs instead of ints because M.innerSize() returns a long. Is this necessary everywhere?

    krylovSchurReturnObject returnObject;

    long k          = k0; // k may need to be adjusted to prevent stagnation
    long N          = M.innerSize();
    long P          = std::min(std::max(2*k,long(20)),N);

    Matrix_cdd V    = Eigen::MatrixXcd::Zero(N,P);
    Matrix_rdd H    = Eigen::MatrixXd::Zero(P,P);
    Vector_cdd v0   = Eigen::VectorXcd::Random(N);
    Vector_cdd v,r;

    Vector_rdd Alpha= Eigen::VectorXd::Zero(P);
    Vector_rdd Beta = Eigen::VectorXd::Zero(P-1);
    Vector_rdd c    = Eigen::VectorXd::Zero(P);
    Vector_rdd d    = Eigen::VectorXd::Zero(P);
    Vector_rdd res  = Eigen::VectorXd::Zero(P);

    Vector_bool isConverged(P);

    Eigen::SelfAdjointEigenSolver<Matrix_rdd> ES;

    long alpha_ind  = 0;
    long d_ind      = 0;
    long n_conv     = 0;
    long sizeV      = 0;
    bool stopAlg    = false;
    bool restarted  = true;
    double normv,normRes=0,alpha=0;

    // Initialise
    v               = M*v0;
    normv           = v.norm();
    v               /= normv;
    normv           = v.norm();
    if (abs(normv-1.0) >= 1e-10){
        v   = v0;
    }

    // Main loop
    for (long iter = 0; iter < maxIter; ++iter) {
        std::cout << iter << std::endl;
        
        // Build Krylov subspace
        for (long i = sizeV; i < P; ++i) {
            std::cout << "\t" << i << std::endl;
            V.col(i)    = v;
            r           = M*v;
            alpha       = v.dot(r).real();

            if (i==1){
                r       = r- alpha*r;
            } else if (restarted) {
                r           = r -projector(V,r,i);
                restarted   = false;
            } else {
                r           = r - alpha*v - normRes*V.col(i-1);
            }

            robustOrthogonalise(V, r, i, normRes, stopAlg);

            if (stopAlg) {
                return returnObject;
            }
            // Store results
            v                   = r;
            Alpha(alpha_ind)    = alpha;
            if (i!=(P-1)) {
                // TODO: Get rid of this branch!
                Beta(alpha_ind) = normRes;
            }
            alpha_ind++;
        }

        // Build H1
        // TODO: How to build H1 when alpha, beta are variably sized? We dont! Just insert straight into H
        H.setZero(); // re-zero


        std::cout << "\t\t" << "A" << std::endl;
        (H.block(d_ind,d_ind,P-d_ind,P-d_ind)).diagonal()       = Alpha.head(alpha_ind);
        std::cout << "\t\t" << "B" << std::endl;
        (H.block(d_ind,d_ind,P-d_ind,P-d_ind)).diagonal(+1)     = Beta.head(alpha_ind-1);
        std::cout << "\t\t" << "C" << std::endl;
        (H.block(d_ind,d_ind,P-d_ind,P-d_ind)).diagonal(-1)     = Beta.head(alpha_ind-1);
        std::cout << "\t\t" << "D" << std::endl;
//        std::cout << H << std::endl;

        std::cout << H.block(0,k,k,1)<< std::endl;
        std::cout << c << std::endl;
                
        if (d_ind!=0) {
            (H.block(0,0,d_ind,d_ind)).diagonal()   = d.head(d_ind);
            H.block(0,k,k,1)                        = c;
            H.block(k,0,1,k)                        = c;
        }
        alpha_ind   = long(0);
        std::cout << "\t\t" << "E" << std::endl;
        // Compute eigenvalues. N.B. Eigen::SelfAdjointEigenSolver sorts eigenvalues in ascending order
        ES.compute(H);
        d               = ES.eigenvalues();

        // Compute residuals
        res             = (normRes*ES.eigenvectors().row(P-1)).cwiseAbs();

        // Count number of converged eigenpairs
        isConverged     = res.head(k).array() < d.head(k).cwiseAbs().array().max(pow(DBL_EPSILON,2.0/3.0))*tol;
        n_conv          = (isConverged).count();

        if (n_conv>=k || iter==maxIter) {
            // We've either reached the required number of converged eigenpairs
            // Or we're out of time
            break;
        } else {
            // Adjust k to prevent stagnation
            k           = k0 + std::min<long>(n_conv,(P-k0)/2); // This should actually be safe because of C++ int division
            if (k==1 && P>3) {
                k       = P/2;
            }
        }

        // Store variables for next iteration
        V.leftCols(k)   = V*ES.eigenvectors().leftCols(k);
        d_ind           = k;
        c               = normRes*ES.eigenvectors().row(0);
        restarted       = true;
        sizeV           = k;
    }

    // Populate returnObject fields
    returnObject.eigenvalues    = d.head(k0);
    returnObject.eigenvectors   = ES.eigenvectors().leftCols(k0);
    returnObject.isConverged    = isConverged.head(k0);

    return returnObject;
}

template <class C>
Matrix_rdd lanczos(const C &M,  const long L) {
    // M is the Hermitian matrix we want to diagonalise
    // M must inherit from Matrix_cdd or Matrix_cds
    // L is the number of Lanczos iterations to perform

    // TODO: Add eigenvectors as return from lanczos dense, e.g. create a struct that contains eigenvalues and eigenvectors. Note that if Tv=\lambda v, then M w = \lambda w, where w = V*v
    // TODO: Check that template works with sparse inputs

    static_assert(std::is_base_of<Matrix_cdd,C>::value || std::is_base_of<Matrix_cds,C>::value,
                  "Input matrix M must inherit from Eigen::MatrixXd or Eigen::MatrixXcd");

    long N = M.innerSize();

    // Matrix to store vector sequence, V
    Matrix_cdd V  = Eigen::MatrixXcd::Zero(N,L);
    // Tridiagonal system to be returned
    Matrix_rdd T = Eigen::MatrixXd::Zero(L,L);

    // Coefficients in tridiagonal matrix T
    Vector_rdd a = Eigen::VectorXd::Zero(L);
    Vector_rdd b = Eigen::VectorXd::Zero(L-1);
    double b0;

    // Random initial vector; should make this deterministic, i.e. fix prng seed
    Vector_cdd v = Eigen::VectorXcd::Random(N);
    // Intermediate vectors
    Vector_cdd w,v_gs;
    // Intermediate values
    double temp_norm;

    // Other flags
    bool is_orthogonal = false;


    b0      = v.norm();
    v       /= b0;

    V.col(0)= v;

    w       = M*v;

    a(0)    = (v.dot(w)).real();
    w       = w-a(0)*v;


    // Lanczos iteration
    for (int i = 1; i < L; ++i) {
        b(i-1)  =  w.norm();

        // Check that norm isn't zero
        if (b(i-1)>DBL_EPSILON){
            // If norm is non-zero, add construct new v
            v = w / b(i-1);
        } else {
            // If norm is zero, use Gram-Schmidt to construct new v orthogonal to existing vector sequence
            while(!is_orthogonal){
                v_gs = Eigen::VectorXcd::Random(N);
                v_gs /= v_gs.norm();
                for (int j = 0; j < i; ++j) {
                    v_gs -= (v.array()*V.array().col(j)).matrix();
                }
                temp_norm = v_gs.norm();
                if (temp_norm>DBL_EPSILON) {
                    v_gs /= temp_norm;
                    is_orthogonal = true;
                }
            }
            v = v_gs;
        }

        V.col(i)    = v;

        w               = M*v;
        a(i)            = (w.dot(v)).real();
        w               = w - a(i)*v - b(i-1)*V.col(i-1);
    }

    // Construct T
    T.diagonal()    = a;
    T.diagonal(+1)  = b;
    T.diagonal(-1)  = b;
    return T;
}


Matrix_rdd lanczos_dense_sparseT(const Matrix_cdd &M,  const long L) {
    // M is the hermitian matrix we want to find eigenvalues for
    // N is the side length of M
    // L is the number of Lanczos iterations to perform
    long N = M.innerSize();

    // Matrix to store vector sequence, V
    Matrix_cdd V  = Eigen::MatrixXcd::Zero(N,L);
    // Tridiagonal system to be diagonalised, T
    Matrix_rds T(L,L);
    T.setZero();
    T.reserve(3*L);

    // Coefficients in tridiagonal matrix T
    //Vector_rdd a = Eigen::VectorXd::Zero(L);
    //Vector_rdd b = Eigen::VectorXd::Zero(L-1);
    double b0;

    // Random initial vector; should make this deterministic, i.e. fix prng seed
    Vector_cdd v = Eigen::VectorXcd::Random(N);
    // Intermediate vectors
    Vector_cdd w,v_gs;
    // Intermediate values
    double temp_norm;

    // Other flags
    bool is_orthogonal = false;


    b0      = v.norm();
    v       /= b0;

    V.col(0)= v;

    w       = M*v;

    T.insert(0,0)   = (v.dot(w)).real();
    w               = w- T.coeff(0,0)*v;


    // Lanczos iteration
    for (int i = 1; i < L; ++i) {
        T.insert(i-1,i)  =  w.norm();

        // Check that norm isn't zero
        if (T.coeff(i-1,i)>DBL_EPSILON){
            // If norm is non-zero, add construct new v
            v = w / T.coeff(i-1,i);
        } else {
            // If norm is zero, use Gram-Schmidt to construct new v orthogonal to existing vector sequence
            while(!is_orthogonal){
                v_gs = Eigen::VectorXcd::Random(N);
                v_gs /= v_gs.norm();
                for (int j = 0; j < i; ++j) {
                    v_gs -= (v.array()*V.array().col(j)).matrix();
                }
                temp_norm = v_gs.norm();
                if (temp_norm>DBL_EPSILON) {
                    v_gs /= temp_norm;
                    is_orthogonal = true;
                }
            }
            v = v_gs;
        }

        V.col(i)    = v;

        w               = M*v;
        T.insert(i,i)   = (w.dot(v)).real();
        w               = w - T.coeff(i,i)*v - T.coeff(i-1,i)*V.col(i-1);
    }

    // Construct T
    T.selfadjointView<Eigen::Lower>() = T.selfadjointView<Eigen::Upper>();
    return Matrix_rdd(T);
}



/*
Matrix_rds lanczos_sparse(const Matrix_cds &M, const long L) {
    long N = M.innerSize();

    // May need to initialise these to zero?
    //Matrix to store vector sequence V
    Matrix_cds V;
    V.reserve(N*L);

    // Tridiagonal system to be diagonalised. There will be at most 3*L-2 non-zero values
    Matrix_rds T;
    T.reserve(3*L-2);

    // Coefficients in tridiagonal matrix T
    */
/*
    Vector_rdd a = Eigen::VectorXd::Zero(L);
    Vector_rdd b = Eigen::VectorXd::Zero(L-1);
    double b0;

    // Random initial vector. We still ought to store this densely
    Vector_cds v = Eigen::VectorXcd::Random(N).sparseView();
    // Intermediate vectors
    Vector_cds w,v_gs;

    // Intermediate values
    double temp_norm;

    // Other flags
    bool is_orthogonal = false;


    b0 = v.norm();
    v /= b0;

    V.col(0)    = v;

    w           = M*v; // assuming all coefficients of M are stored

    //a(0) = (v.dot(w)).real();
    T.insert(0,0)   = (v.dot(w)).real();
    w               = w - a(0)*v;

    // Lanczos iteration
    for (int i = 1; i < L; ++i) {
        //b(i-1)  = w.norm();
        T.insert(i-1,i)   = w.norm();

        // Check that norm isn't zero
        if (b(i-1)>DBL_EPSILON){
            // If norm is non-zero, add construct new v
            v = w / T.coeff(i-1,i);
        } else {
            // If norm is zero, use Gram-Schmidt to construct new v orthogonal to existing vector sequence
            while(is_orthogonal==false){
                v_gs = Eigen::VectorXcd::Random(N).sparseView();
                v_gs /= v_gs.norm();
                for (int j = 0; j < i; ++j) {
                    v_gs -= (v.cwiseProduct(V.col(j)));
                }
                temp_norm = v_gs.norm();
                if (temp_norm>DBL_EPSILON) {
                    v_gs /= temp_norm;
                    is_orthogonal = true;
                }
            }
            v = v_gs;
        }

        V.col(i)    = v;

        w           = M*v;
        //a(i)        = (v.dot(w)).real();
        T.insert(i,i)       = (v.dot(w)).real();
        //w                   = w - a(i)*v - b(i-1)*V.col(i-1);
        w                   = w - T.coeff(i,i)*v - T.coeff(i-1,i)*V.col(i-1);
    }
    *//*

    // Construct T
    T.selfadjointView<Eigen::Lower>() = T.selfadjointView<Eigen::Upper>();

    return T;

}
*/

Vector_rdd eigs(const Matrix_cdd &M, const int L) {
    Matrix_rdd T = lanczos(M, L);

    Eigen::SelfAdjointEigenSolver<Matrix_rdd> eigensolver(T);
    //std::cout << eigensolver.eigenvectors().cwiseAbs() << std::endl;
    return eigensolver.eigenvalues().real();
}

/*
Vector_rdd eigs_sparse(const Matrix_cds &M, const int L) {
    Matrix_rdd T = lanczos_sparse(M, L);

    Eigen::EigenSolver<Matrix_rds> eigensolver(T);
    return eigensolver.eigenvalues().real();
}
*/

Matrix_cdd gen_hermitian_dense(const int n)
{
	const Matrix_cdd mat = Eigen::MatrixXcd::Random(n, n);
	return mat + mat.adjoint();
}

/*
Matrix_cds gen_hermitian_sparse(const int n)
{
    const Matrix_cds mat = Eigen::MatrixXcd::Random(n, n).sparseView();
    return mat + mat.adjoint();
}
*/


Matrix_cdd gen_dense_mat()
{
	Matrix_cdd mat(3,3);
	mat(0,0) 	= 1.0 + 1i;
	mat(0,1) 	= 1.0 - 1i;
	mat(0,2) 	= -1.0 + 1i;
	mat(1,0) 	= -1.0 - 1i;
	mat(1,1) 	= 2.0 + 7.0*1i;
	mat(1,2) 	= 3.0 + 8.0*1i;
	mat(2,0) 	= 4.0 + 9.0*1i;
	mat(2,1) 	= 5.0 + 10.0*1i;
	mat(2,2) 	= 6.0 + 11.0*1i;
	return mat;
}



int main()
{
	int N = 10, L=6;

	// DENSE
	const Matrix_cdd M0 = gen_hermitian_dense(N);

	// Explicit
    Vector_rdd evals_exp;
    double eval_exp_min;
    auto start = std::chrono::high_resolution_clock::now();
	Eigen::ComplexEigenSolver<Matrix_cdd> eigensolver(M0);
    evals_exp = eigensolver.eigenvalues().real();
    eval_exp_min = evals_exp.minCoeff();
    auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish-start;

    // Implicit
    Vector_rdd evals_lanczos;
    double eval_lanczos_min;
    auto start0 = std::chrono::high_resolution_clock::now();
    evals_lanczos = eigs(M0, L);
    eval_lanczos_min = evals_lanczos.minCoeff();
    auto finish0 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed0 = finish0-start0;


    // krylovSchur(const C &M, const long k0, const long maxIter = 300, const double tol = 1.0e-14) {
    
    // Implicit: Krylov-Schur

    long k0 = 6;
    double eval_krylovschur_min;
    auto start1 = std::chrono::high_resolution_clock::now();
    krylovSchurReturnObject ksro =  krylovSchur(M0, 6);
    eval_krylovschur_min = ksro.eigenvalues.minCoeff();
    auto finish1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = finish1-start1;

    std::cout << "DENSE" << std::endl;
    std::cout << "Exact:\t" << eval_exp_min << "\tTime taken:\t" << elapsed.count() << std::endl;
    std::cout << "Lanczos:\t" << eval_lanczos_min << "\tTime taken:\t"<< elapsed0.count()  << std::endl;
    std::cout << "Diff:\t" << abs(eval_lanczos_min-eval_exp_min) << std::endl;
    std::cout << "Krylov-Schur:\t" << eval_krylovschur_min << "\tTime taken:\t" << elapsed1.count() << std::endl;
    std::cout << "Diff:\t" << abs(eval_krylovschur_min-eval_exp_min) << std::endl;



    // SPARSE
    //const Matrix_cds M0_sparse = gen_hermitian_sparse(N);

/*    // Explicit
    Vector_rdd evals_exp_sparse;
    double eval_exp_min_sparse;
    Eigen::ComplexEigenSolver<Matrix_cds> eigensolver_sparse(M0_sparse);
    evals_exp = eigensolver.eigenvalues().real();
    eval_exp_min = evals_exp.minCoeff();

    // Implicit
    Vector_rdd evals_lanczos_sparse;
    double eval_lanczos_min_sparse;
    evals_lanczos_sparse = eigs_sparse(M0_sparse, L);
    eval_lanczos_min_sparse = evals_lanczos.minCoeff();


    std::cout << "SPARSE" << std::endl;
    std::cout << "Exact:\t" << eval_exp_min_sparse << std::endl;
    std::cout << "Lanczos:\t" << eval_lanczos_min_sparse << std::endl;
    std::cout << "Diff:\t" << abs(eval_lanczos_min_sparse-eval_exp_min_sparse) << std::endl;*/
/*
    Vector_cdd V = Eigen::VectorXcd::Random(4);
    Vector_cdd V2 = Eigen::VectorXcd::Random(3);
    Matrix_cdd M = Eigen::MatrixXcd::Random(4,4);
    std::cout << std::endl << V << std::endl << std::endl;
    V /= V.norm();
    std::cout << V << std::endl << std::endl;
    std::cout << M << std::endl;
    M.diagonal() = V;
    M.diagonal(+1) = V2;
    std::cout << M << std::endl;
*/
}


