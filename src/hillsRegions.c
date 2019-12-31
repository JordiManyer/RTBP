
#include "../myHeader.h"

void hillRegion(HR reg , double C , double mu){
    double x1,x2,x3;
    double C1,C2,C3,C4;
    double L[4];

    // We first look for L1,2,3,4
    x1 = solveL1(mu , 1e-14);
    x2 = solveL2(mu , 1e-14);
    x3 = solveL3(mu , 1e-14);
    for (int i = 0 ; i < 4 ; ++i) L[i] = 0;
    L[0] = x1; C1 = Jacobi(L , mu);
    L[0] = x2; C2 = Jacobi(L , mu);
    L[0] = x3; C3 = Jacobi(L , mu);
    L[0] = mu -1.0/2.0; L[1] = sqrt(3.)/2.0; C4 = Jacobi(L , mu);

    //


}

void FHill(double F, double Fx , double Fy , double* x , double mu , double C){
    double r1, r2, Omega, Ox, Oy;
    r1 = sqrt(pow(x[0]-mu , 2) + pow(x[1],2));
    r2 = sqrt(pow(x[0]-mu+1 , 2) + pow(x[1],2));

    Omega = 0.5*(x[0]*x[0] + x[1]*x[1]) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu);
    Ox = x[0] - ((1.0-mu)*(x[0]-mu))/(pow(r1,3)) - (mu*(x[0]-mu+1))/(pow(r2,3));
    Oy = x[1]*( 1 - (1-mu)/(pow(r1,3)) - mu/(pow(r2,3)) );

    F = 2.0 * Omega - C;
    Fx = 2.0 * Ox;
    Fy = 2.0 * Oy;
}