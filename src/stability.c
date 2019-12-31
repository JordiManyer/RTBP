

#include "../myHeader.h"


void getStabilityFlows(MMat M){
    int flag, dim;
    gsl_matrix* auxMat;
    gsl_eigen_nonsymmv_workspace * w;
    gsl_complex lambda; double realpart;


    // Init workspace
    dim = M.n;
    w = gsl_eigen_nonsymmv_alloc (dim);
    gsl_eigen_nonsymmv_params(1 , w);
    auxMat = gsl_matrix_alloc (dim, dim);

    // Copy matrix into it'c copy (needed since the solver destroys the matrix)
    flag = gsl_matrix_memcpy(auxMat, M.M);

    // Calculate eigenvectors and eigenvalues
    flag = gsl_eigen_nonsymmv(auxMat, M.eval, M.evec, w);
    gsl_eigen_nonsymmv_free (w);

    // Sort eigenvalues
    for (int i = 0; i < dim; i++){
        lambda = gsl_vector_complex_get (M.eval, i);
        realpart = GSL_REAL(lambda);
        if (realpart < 0.0) M.stable[i] = 1;
        else if (realpart > 0.0) M.stable[i] = 2;
        else if (fabs(realpart) < 1.e-10) M.stable[i] = 3;
        else M.stable[i] = -1;
    }
}

void getStabilityMaps(MMat M){
    int flag, dim;
    gsl_matrix* auxMat;
    gsl_eigen_nonsymmv_workspace * w;
    gsl_complex lambda; double realpart;

    // Init workspace
    dim = M.n;
    w = gsl_eigen_nonsymmv_alloc (dim);
    gsl_eigen_nonsymmv_params(1 , w);
    auxMat = gsl_matrix_alloc (dim, dim);

    // Copy matrix into it'c copy (needed since the solver destroys the matrix)
    flag = gsl_matrix_memcpy(auxMat, M.M);

    // Calculate eigenvectors and eigenvalues
    flag = gsl_eigen_nonsymmv(auxMat, M.eval, M.evec, w);
    gsl_eigen_nonsymmv_free (w);

    // Sort eigenvalues
    for (int i = 0; i < dim; i++){
        lambda = gsl_vector_complex_get (M.eval, i);
        realpart = gsl_complex_abs(lambda);
        if (realpart < 1.0) M.stable[i] = 1;
        else if (realpart > 1.0) M.stable[i] = 2;
        else if (fabs(realpart - 1.0) < 1.e-10) M.stable[i] = 3;
        else M.stable[i] = -1;
    }
}

void getParametrization(gsl_matrix* cj ,int nPoints,int dim, double T, double* x0, MMat M , double s) {
    int iStable, iUnst;
    double eigStable, eigUnst;
    gsl_complex lambda; double realpart;
    gsl_vector* v0; double aux;

    // Get stable and unstable eigenvalues and their positions
    iStable = -1; iUnst = -1;
    eigStable = 100.0; eigUnst = 0.0;
    for (int i = 0; i < dim ; ++i) {
        lambda = gsl_vector_complex_get (M.eval, i);
        if (GSL_IMAG(lambda) > 1.e-10) printf("WARNING: Eigenvalue is not full real!!");
        realpart = GSL_REAL(lambda);

        if (fabs(realpart) < fabs(eigStable)) {iStable = i; eigStable = realpart;}
        if (fabs(realpart) > fabs(eigUnst)) {iUnst = i; eigUnst = realpart;}
    }

    // UNSTABLE MANIFOLD: Initial values saved in columns 0,...,2*nPoints-1 of cj
    v0 = gsl_vector_calloc(dim);
    for (int i = 0; i < dim; ++i) {
        if (GSL_IMAG(gsl_matrix_complex_get(M.evec , iUnst , i)) > 1.e-10) printf("WARNING: Eigenvector is not full real!!");
        gsl_vector_set(v0 , i , GSL_REAL(gsl_matrix_complex_get(M.evec , iUnst , i)));
    }

    for (int j = 0; j < nPoints; ++j){
        aux = s * pow(eigUnst , -(double)(j)/(double)(nPoints-1));
        for (int i = 0; i < dim; ++i) {
            gsl_matrix_set(cj , i , j , x0[i] + aux * gsl_vector_get(v0 , i));
            gsl_matrix_set(cj , i , nPoints + j , x0[i] - aux * gsl_vector_get(v0 , i));
        }
    }
    gsl_vector_free(v0);

    // STABLE MANIFOLD: Initial values saved in columns 2*nPoints,...,4*nPoints-1 of cj
    v0 = gsl_vector_calloc(dim);
    for (int i = 0; i < dim; ++i) {
        if (GSL_IMAG(gsl_matrix_complex_get(M.evec , iStable , i)) > 1.e-10) printf("WARNING: Eigenvector is not full real!!");
        gsl_vector_set(v0 , i , GSL_REAL(gsl_matrix_complex_get(M.evec , iStable , i)));
    }

    for (int j = 0; j < nPoints; ++j){
        aux = s * pow(eigStable , -(double)(j)/(double)(nPoints-1));
        for (int i = 0; i < dim; ++i) {
            gsl_matrix_set(cj , i , 2*nPoints + j , x0[i] + aux * gsl_vector_get(v0 , i));
            gsl_matrix_set(cj , i , 3*nPoints + j , x0[i] - aux * gsl_vector_get(v0 , i));
        }
    }
    gsl_vector_free(v0);
}


void arrayToMatrix(gsl_matrix* M , Orbit_Var orb){
    int nPoints = orb.nt;

    gsl_matrix_set(M, 0 , 0 , orb.x5[nPoints-1]);
    gsl_matrix_set(M, 0 , 1 , orb.x6[nPoints-1]);
    gsl_matrix_set(M, 0 , 2 , orb.x7[nPoints-1]);
    gsl_matrix_set(M, 0 , 3 , orb.x8[nPoints-1]);

    gsl_matrix_set(M, 1 , 0 , orb.x9[nPoints-1]);
    gsl_matrix_set(M, 1 , 1 , orb.x10[nPoints-1]);
    gsl_matrix_set(M, 1 , 2 , orb.x11[nPoints-1]);
    gsl_matrix_set(M, 1 , 3 , orb.x12[nPoints-1]);

    gsl_matrix_set(M, 2 , 0 , orb.x13[nPoints-1]);
    gsl_matrix_set(M, 2 , 1 , orb.x14[nPoints-1]);
    gsl_matrix_set(M, 2 , 2 , orb.x15[nPoints-1]);
    gsl_matrix_set(M, 2 , 3 , orb.x16[nPoints-1]);

    gsl_matrix_set(M, 3 , 0 , orb.x17[nPoints-1]);
    gsl_matrix_set(M, 3 , 1 , orb.x18[nPoints-1]);
    gsl_matrix_set(M, 3 , 2 , orb.x19[nPoints-1]);
    gsl_matrix_set(M, 3 , 3 , orb.x20[nPoints-1]);
}

