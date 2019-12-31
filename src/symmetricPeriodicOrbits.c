
#include "../myHeader.h"

double computeYPrime(double C , double x0);
double computeXprime(double C , double x0);
double mybissection( double (*func)(double , double) , double a , double b , double C , double tol);


double findSymmetricPO(double C , double x0, double xstep) {
    int k;
    double x1, xp1, x2, xp2;
    double xfinal;
    double (*optfunc)(double , double);

    // Initial points
    x1 = x0 + xstep / 2.0;
    xp1 = computeXprime(C, x1);
    x2 = x0 + xstep;
    xp2 = computeXprime(C, x2);

    if (verbose) {
        printf("Computations for C = %g \n" , C);
        printf("Starting points: %g %g \n" , x1 , x2);
        printf("xprime values = %g %g \n" , xp1 , xp2);
    }

    k = 0;
    while (xp1 * xp2 > 0){
        x1 = x2;
        xp1 = xp2;

        x2 = x2 + xstep;
        xp2 = computeXprime(C, x2);
        k = k + 1;

        if (verbose) {
            printf("Iteration %d \n" , k);
            printf("New iterate: %g \n" , x2);
            printf("New xprime value = %g \n" , xp2);
        }
    }


    if(verbose) printf("Initial condition found: %g %g \n" , x1 , x2);
    optfunc = computeXprime;
    xfinal = mybissection(optfunc, x1, x2, C, 1e-15);
    return xfinal;
}


void computeSymmetricPO(Orbit orb , double C , double x0 , int nPoints){
    Point p;
    double (*event)(double*);
    double* ic; double x[4];

    event = yequalzero;
    x[0] = x0; x[1] = 0;
    x[2] = 0; x[3] = computeYPrime(C,x0);
    ic = x;

    // Compute final time
    p = computePoincare(ic , 0.0 , 1 , event , 2);
    // Compute orbit up to final time
    computeOrbitNpoints(orb , ic , 0.0 , p.t , nPoints, 1);
}


double computeYPrime(double C , double x0){
    extern MY_FLOAT mu;
    double r1 , r2 , Omega , yp0;

    r1 = (x0-mu);
    r2 = (x0-mu+1);
    Omega = 0.5*(x0*x0) + (1-mu)/r1 + mu/r2 + 0.5*mu*(1-mu);

    yp0 = -sqrt(2.0*Omega - C);
    return yp0;
}


double computeXprime(double C , double x0){
    Point p;
    double (*event)(double*);
    double* ic; double x[4];

    event = yequalzero;
    x[0] = x0; x[1] = 0;
    x[2] = 0; x[3] = computeYPrime(C,x0);
    ic = x;

    p = computePoincare(ic , 0.0 , 1 , event , 1);
    return p.x[2];
}


double mybissection( double (*func)(double , double) , double a , double b , double C , double tol) {
    double mid, fa, fb, fmid;
    int done, k;

    fa = func(C, a);
    fb = func(C, b);
    if (fa * fb > 0) {
        return -1;
    }
    if (verbose) {
        printf("Starting Bissection: \n");
        printf("a: %g , b: %g \n", a , b);
        printf("fa: %g , fb: %g \n", fa , fb);
    }

    done = 0; k = 0;
    mid = (a+b)/2.0;
    while (! done) {
        mid = (a + b) / 2.0;
        fmid = func(C, mid);

        if (fmid * fa < 0) {
            b = mid;
            fb = fmid;
        } else {
            a = mid;
            fa = fmid;
        }

        if (verbose) {
            printf("Iteration %d: \n", k);
            printf("a: %g , b: %g , mid: %g \n", a , b , mid);
            printf("fa: %g , fb: %g , fmid: %g \n", fa , fb, fmid);
            printf("abs(a-b): %g \n", fabs(a-b));
        }

        if (fabs(a-b) < tol) done = 1;
        ++k;
    }
    if (verbose) printf("Solution found: %g \n", mid);


    return mid;
}
