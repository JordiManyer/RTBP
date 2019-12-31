#include "../myHeader.h"

double F1(double x, MY_FLOAT mu);
double F2(double x, MY_FLOAT mu);
double F3(double x, MY_FLOAT mu);

// SOLVERS
MY_FLOAT solveL1(MY_FLOAT mu , double tol)
{
    double xkm1 = pow(mu/(3.*(1.-mu)) , 1./3.);
    double xk = F1(xkm1,mu);
    while (fabs(xkm1 -xk)  > tol)
    {
        xkm1 = xk;
        xk = F1(xk,mu);
    }
    return mu - 1 - xk;
}

MY_FLOAT solveL2(MY_FLOAT mu , double tol)
{
    double xkm1 = pow(mu/(3.*(1.-mu)) , 1./3.);
    double xk = F2(xkm1,mu);
    while (fabs(xkm1 -xk)  > tol)
    {
        xkm1 = xk;
        xk = F2(xk,mu);
    }
    return mu - 1 + xk;
}

MY_FLOAT solveL3(MY_FLOAT mu , double tol)
{
    double xkm1 = 1.0 - 7.0/12.0 * mu;
    double xk = F3(xkm1,mu);
    while (fabs(xkm1 -xk)  > tol)
    {
        xkm1 = xk;
        xk = F3(xk,mu);
    }
    return mu + xk;
}


// FIXED POINT FUNCTIONS
double F1(double x, MY_FLOAT mu)
{
    double y;
    y = pow( mu * pow(1 - x ,2) / (3 - 2*mu - x*(3 - mu - x)) , 1.0/3.0 );
    return y;
}

double F2(double x, MY_FLOAT mu)
{
    double y;
    y = pow( mu * pow(1 + x ,2) / (3 - 2*mu + x*(3 - mu + x)) , 1.0/3.0 );
    return y;
}

double F3(double x, MY_FLOAT mu)
{
    double y;
    y = pow( (1.0 - mu) * pow(1.0 + x , 2.0) / (1.0 + 2.0*mu + x*(2.0 + mu + x)) , 1.0/3.0 );
    return y;
}





