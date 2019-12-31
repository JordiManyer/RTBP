#include "../myHeader.h"

// Regular equations
void RTBP_F(double* res , double* x,double mu){
    double r1, r2, Ox, Oy;

    r1 = sqrt(pow(x[0]-mu , 2) + pow(x[1],2));
    r2 = sqrt(pow(x[0]-mu+1 , 2) + pow(x[1],2));
    Ox = x[0] - ((1.0-mu)*(x[0]-mu))/(pow(r1,3)) - (mu*(x[0]-mu+1))/(pow(r2,3));
    Oy = x[1]*( 1 - (1-mu)/(pow(r1,3)) - mu/(pow(r2,3)) );

    res[0] = x[2]; res[1] = x[3]; res[2] = 2.0*x[3] + Ox;  res[3] = -2.0*x[2] + Oy;
}


// Variational Matrix
void RTBP_VarMat(gsl_matrix* A, double* x,double mu){
    double r1, r2, Oxx, Oxy, Oyy;

    r1 = sqrt(pow(x[0]-mu , 2) + pow(x[1],2));
    r2 = sqrt(pow(x[0]-mu+1 , 2) + pow(x[1],2));

    Oxx = 1 + ((3./4.)*(1.-mu)*pow(2.*x[0]-2.*mu,2))/pow(r1,5) - (1.-mu)/pow(r1,3) + ((3./4.)*(mu)*pow(2.*x[0]-2.*mu + 2.,2))/pow(r2,5) - mu/pow(r2,3);
    Oyy = 1 + ((3./1.)*(1.-mu)*pow(x[1],2))/pow(r1,5) - (1-mu)/pow(r1,3) + ((3./1.)*(mu)*pow(x[1],2))/pow(r2,5) - mu/pow(r2,3);
    Oxy = (3./2.)*(1.-mu)*(2.*x[0]-2.*mu)*x[1]/pow(r1,5) + (3./2.)*(1.-mu)*(2.*x[0]-2.*mu+2.)*x[1]/pow(r2,5);

    A = gsl_matrix_calloc (4, 4); // Init to zero

    // Equations
    gsl_matrix_set(A, 0, 2, 1.0);
    gsl_matrix_set(A, 1, 3, 1.0);

    gsl_matrix_set(A, 2, 3, 2.0);
    gsl_matrix_set(A, 2, 0, Oxx);
    gsl_matrix_set(A, 2, 1, Oxy);

    gsl_matrix_set(A, 3, 2, -2.0);
    gsl_matrix_set(A, 3, 0, Oxy);
    gsl_matrix_set(A, 3, 2, Oyy);
}

// Jacobi first integral
double Jacobi(double* x, double mu){
    double r1, r2, Omega, C;

    r1 = sqrt(pow(x[0]-mu , 2) + pow(x[1],2));
    r2 = sqrt(pow(x[0]-mu+1 , 2) + pow(x[1],2));

    Omega = 0.5*(x[0]*x[0] + x[1]*x[1]) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu);

    C = 2.0*Omega - (x[2]*x[2] + x[3]*x[3]);
    return C;
}
