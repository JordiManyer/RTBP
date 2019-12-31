
#include "myHeader.h"

const int orbitMaxLength = 1000; // Maximum number of pts in in a single orbit

// Verbose Levels:
// 2 -> debugging + program steps + results
// 1 -> program steps + results ,
// 0 -> only results
const int verbose = 0;

double PoincareValue;
MY_FLOAT mu;                     // Global variable mu

// DIFFERENT PREPARED ROUTINES:
void demo();
void Assignment11();
void Assignment12();
void Assignment13();

// MAIN PROGRAM
int main(int argc, char **argv)
{
    //demo();
    Assignment12();
    exit(0);
}


// --------------------------------------------------------
// DEMO ROUTINE TO TEST FUNCTIONALITIES
void demo(){
    // L3 vars
    double    L3tol = 1.e-14;
    double    L3[4];

    // Integration params
    double* x0;
    double tstart, tend, tstep;
    int dir;

    // Poincare orbit
    int nCrossings;
    double (*pcfunc)(double*);

    // Structs
    Orbit orb;
    Point p;

    // Calculate L3
    mu = 0.1;
    L3[0] = solveL3(mu , L3tol);
    L3[1] = L3[2] = L3[3] = 0.0;


    // Integrate normal orbit near L3
    tstart = 0.0;
    tend = 1.0;
    tstep = 0.1;
    dir = 1;
    x0 = L3; x0[0] = x0[0] + 0.1;

    orb = computeOrbit(x0 , tstart , tend , tstep , dir);
    printOrbit(orb);

    // Poincare Section up to first crossing with x == 0
    pcfunc = xequalzero;
    nCrossings = 1;
    p = computePoincare(x0,tstart,dir,pcfunc,nCrossings);
    printPoint(p);

}

// --------------------------------------------------------
// ASSIGNEMENT 11
void Assignment11(){
    double C, x0, xPO, xstep;
    double Cstart, Cend, Cresolution;
    double L3tol, xL3;
    int numPO, nPoints;
    FILE *fptr;
    Orbit orb;

    // Allocate space for the periodic orbit
    nPoints = 200;
    double t[nPoints] , x[nPoints] , y[nPoints], xp[nPoints], yp[nPoints];
    orb.nt = nPoints; orb.t = t; orb.x1 = x; orb.x2 = y; orb.x3 = xp; orb.x4 = yp;

    // Set range for C
    Cstart = 3.189; Cend = 2.1; Cresolution = 1.e-3;

    // Open the output file
    fptr = fopen("Ass11.out","w");

    // Calculate L3
    mu = 0.1; L3tol = 1.e-15;
    xL3 = solveL3(mu , L3tol);

    // Find periodic orbits:
    numPO = ceil((Cstart-Cend)/Cresolution); fprintf(fptr , "%d \n" , numPO);
    xPO = xL3;
    for (int i = 0 ; i <= numPO ; ++i) {
        if (i == 0) xstep = 1.e-5;
        else xstep = 1.e-4;

        x0 = xPO;
        C = Cstart - (Cstart-Cend)*i/numPO;
        xPO = findSymmetricPO(C , x0, xstep);
        printf("PO orbit found: C=%g , x=%g \n", C , xPO);

        // Compute the orbit and print it to file
        computeSymmetricPO(orb , C , xPO , nPoints);
        fprintf(fptr,"%g\n",C); // Print C
        tofileOrbit(fptr , orb); // Print orbit
    }
    fclose(fptr);
}

// --------------------------------------------------------
// ASSIGNEMENT 12
void Assignment12(){
    double T, Thalf, x0[20], stable[4];
    double s, x0Manifold[4];
    double (*event)(double*);
    int dim, nPoints , nPointsPrint;
    FILE *fptr;
    Orbit_Var orb; Orbit orb2;
    MMat M;
    gsl_matrix* cj;
    gsl_matrix* vtheta;

    // Initial conditions for the PO
    mu = 0.01;
    Thalf = 3.114802556760205;
    x0[0] = 1.033366313746765;
    x0[1] = 0.0; x0[2] = 0.0;
    x0[3] = -0.05849376854515592;
    T = 2.0 * Thalf;

    // Initial conditions for the variational matrix (Id)
    for (int i = 4; i < 20; ++i) x0[i] = 0;
    x0[4] = 1.0; x0[9] = 1.0; x0[14] = 1.0; x0[19] = 1.0;

    // Allocate space for the periodic orbit
    nPoints = 30;
    double t[nPoints] , x[nPoints] , y[nPoints], xp[nPoints], yp[nPoints];
    double x5[nPoints] , x6[nPoints], x7[nPoints], x8[nPoints];
    double x9[nPoints] , x10[nPoints], x11[nPoints], x12[nPoints];
    double x13[nPoints] , x14[nPoints], x15[nPoints], x16[nPoints];
    double x17[nPoints] , x18[nPoints], x19[nPoints], x20[nPoints];

    orb.nt = nPoints; orb.t = t; orb.x1 = x; orb.x2 = y; orb.x3 = xp; orb.x4 = yp;
    orb.x5 = x5; orb.x6 = x6; orb.x7 = x7; orb.x8 = x8;
    orb.x9 = x9; orb.x10 = x10; orb.x11 = x11; orb.x12 = x12;
    orb.x13 = x13; orb.x14 = x14; orb.x15 = x15; orb.x16 = x16;
    orb.x17 = x17; orb.x18 = x18; orb.x19 = x19; orb.x20 = x20;

    // Calculate the orbit with first variational equations
    computeVariationalsNpoints(orb, x0 , 0.0 , T , nPoints , 1);
    double errnorm = sqrt( pow(x0[0]-orb.x1[nPoints-1],2) + pow(x0[1]-orb.x2[nPoints-1],2) + pow(x0[2]-orb.x3[nPoints-1],2) + pow(x0[3]-orb.x4[nPoints-1],2));
    double xfinal[4]; xfinal[0] = orb.x1[nPoints-1]; xfinal[1] = orb.x2[nPoints-1]; xfinal[2] = orb.x3[nPoints-1]; xfinal[3] = orb.x4[nPoints-1];
    printf("Error in the periodic orbit: %g \n" , errnorm);
    printf("Error in the final time: %6.12g\n" , fabs(T-orb.t[nPoints-1]));
    printf("Error in the Jacobi invariant : %6.12g\n" , fabs(Jacobi(x0 , mu) - Jacobi(xfinal,mu)));

    // Create and get the stability of the monodromy matrix of the PO
    dim = 4;
    M.n = dim;
    M.M = gsl_matrix_calloc(dim, dim);
    M.eval = gsl_vector_complex_calloc (dim);
    M.evec = gsl_matrix_complex_calloc (dim, dim);
    M.stable = stable;

    arrayToMatrix(M.M , orb);
    getStabilityMaps(M);

    // Safety Check: Calculate determinant of the monodromy matrix (should be 1.0)
    gsl_complex lambda, detM;
    GSL_SET_COMPLEX(&detM , 1.0 , 0.0);
    printf("Eigenvalues of Monodromy matrix , stability (1-> stable, 2-> unstable):\n");
    for (int i = 0; i < dim; ++i){
        lambda = gsl_vector_complex_get(M.eval , i);
        printf("%g + i %g, %d\n" , GSL_REAL(lambda) , GSL_IMAG(lambda) , M.stable[i]);
        detM = gsl_complex_mul(detM , lambda);
    }
    printf("\n");
    printf("Determinant of Monodromy matrix: %g + i %g\n" , GSL_REAL(detM) , GSL_IMAG(detM));

    // Parametrization of the manifold
    s = 1.e-4; // How far are we going from the PO
    cj = gsl_matrix_calloc(dim, 4*nPoints); // Initial points for the 1D manifolds

    /* After getParametrization, cj contains
     * - Initial values for the unstable manifold saved in columns 0,...,2*nPoints-1
     * - Initial values for the stable manifold saved in columns 2*nPoints,...,4*nPoints-1
     */
    getParametrization(cj , nPoints, dim, T, x0, M , s);

    // Open the output file
    fptr = fopen("Ass12.out","w");
    fprintf(fptr , "%d \n" , nPoints);

    // Allocate space for the manifolds
    nPointsPrint = 1000; // number of points each orbit will have
    double t2[nPointsPrint] , x2[nPointsPrint] , y2[nPointsPrint], xp2[nPointsPrint], yp2[nPointsPrint];
    orb2.nt = nPointsPrint; orb2.t = t2; orb2.x1 = x2; orb2.x2 = y2; orb2.x3 = xp2; orb2.x4 = yp2;

    // Integrate unstable manifold initial conditions forward
    printf("Integrating Unstable manifolds.\n");
    PoincareValue = 0.2; event = xequalvalue;
    for (int j = 0; j < 2*nPoints ; ++j) {
        for (int i = 0; i < dim; ++i) {x0Manifold[i] = gsl_matrix_get(cj , i , j);} // get initial cond
        Point p;
        p = computePoincare(x0Manifold , 0.0 , 1 , event , 1); // compute final time
        computeOrbitNpoints(orb2, x0Manifold , 0.0 , p.t , nPointsPrint , 1); // compute manifold
        tofileOrbit(fptr , orb2);
    }

    // Integrate stable manifold initial conditions backwards
    printf("Integrating stable manifolds.\n");
    PoincareValue = 0.2; event = xequalvalue;
    for (int j = 2*nPoints; j < 4*nPoints ; ++j) {
        for (int i = 0; i < dim; ++i) {x0Manifold[i] = gsl_matrix_get(cj , i , j);} // get initial cond
        Point p;
        p = computePoincare(x0Manifold , 0.0 , -1 , event , 1); // compute final time
        computeOrbitNpoints(orb2, x0Manifold , 0.0 , p.t , nPointsPrint , -1); // compute manifold
        tofileOrbit(fptr , orb2);
    }
    printf("Done with asignment 12!!!.\n");
}

void Assignment13(){

}