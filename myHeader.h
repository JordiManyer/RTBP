
#ifndef _MYHEADER_H
#define _MYHEADER_H


// EXTERNAL LIBRARIES:
#include "src/taylor.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

// GLOBAL PARAMETERS:
extern const int verbose; // Different levels of debug output
extern const int orbitMaxLength; // Maximum lengths of the orbits
extern double PoincareValue; // Parameter for the PC section functions
extern double PCTol; // Tolerance for the PC timestep

/*****************************************************************************/
/*****************************************************************************/

// STRUCTURES:
struct orbit{ // Regular orbit
    int nt;
    double *t;
    double *x1;
    double *x2;
    double *x3;
    double *x4;
};
typedef struct orbit Orbit;

/*****************************************************************************/

struct point{ // Single 1+4 dimensional point
    double t;
    double x[4];
};
typedef struct point Point;

/*****************************************************************************/

struct monodromyMatrix{ // Monodromy matrix + stability features
    int n; // dimension of the matrix.
    gsl_matrix * M; // Monodromy Matrix
    gsl_vector_complex * eval; // Vector of eigenvalues
    gsl_matrix_complex * evec; // Matrix which columns are eigenvectors
    int* stable; // Classification vector: 1->stable eigenvalue, 2->Unstable eigenvalue, 3->Pure complex eigenvalue
};
typedef struct monodromyMatrix MMat;

/*****************************************************************************/

struct orbit_var{ // Orbit with variational variables
    int nt;
    double *t;
    // Regular vars: x,y,x',y'
    double *x1; double *x2; double *x3; double *x4;
    // Variational vars: dxi/dxj
    double *x5; double *x6; double *x7; double *x8;
    double *x9; double *x10; double *x11; double *x12;
    double *x13; double *x14; double *x15; double *x16;
    double *x17; double *x18; double *x19; double *x20;
};
typedef struct orbit_var Orbit_Var;

/*****************************************************************************/

struct hill_region{
    int n;
    double* x;
    double* y;
};
typedef struct hill_region HR;

struct hill_region_vector{
    double C;
    int m;
    HR* regions;
};
typedef struct hill_region_vector HRVec;

/*****************************************************************************/
/*****************************************************************************/

// FUNCTION DECLARATIONS BY FILE:
// equPoints.c
MY_FLOAT solveL1(MY_FLOAT mu , double tol); // finds the equilibrium point L1
MY_FLOAT solveL2(MY_FLOAT mu , double tol); // finds the equilibrium point L2
MY_FLOAT solveL3(MY_FLOAT mu , double tol); // finds the equilibrium point L3

// orbitIntegration.c
Orbit computeOrbit(double* x0 , double tstart , double tend , double tstep , int dir); // computes an Orbit of the rtbp
void computeOrbitNpoints(Orbit orb, double* x0 , double tstart , double tend , int nPoints , int dir); // computes an Orbit of the rtbp ensuring nPoints orbit points.
void computeVariationalsNpoints(Orbit_Var orb, double* x0 , double tstart , double tend , int nPoints , int dir);


// poincare.c
Point computePoincare(double* x0 , double tstart , int dir , double (*event)(double*), int nCrossings); // computes the poincaré map of the rtbp
void computeOrbitPoincare(Orbit orb, double* x0 , double tstart , int nPoints , int dir , double (*event)(double*)); // computes the Poincaré map orbit


// poincareSections.c, stopping conditions for computePoincare
double xequalzero(double *x);      // Stop condition x==0
double yequalzero(double *x);      // Stop condition y==0
double xequalvalue(double *x);     // Stop condition x==PoincareValue
double xprimeequalzero(double *x); // Stop condition x'==0

// symmetricPeriodicOrbits.c
double findSymmetricPO(double C , double x0, double xstep); // Finds symmetric periodic orbits around an equilibrium point.
void computeSymmetricPO(Orbit orb , double C , double x0 , int nPoints); //Returns the symmetric orbit found by findSymmetricPO()

// rtbpEqs.c
void RTBP_F(double* res, double* x, double mu);        // Calculate the rtbp equations F(x)
void RTBP_VarMat(gsl_matrix* A, double* x, double mu); // Calculates variational matrix
double Jacobi(double* x, double mu);                   // Calculates the Jacobi first integral

// Stability.c
void getStabilityFlows(MMat M);
void getStabilityMaps(MMat M);
void getParametrization(gsl_matrix* cj ,int nPoints,int dim, double T, double* x0, MMat M , double s);
void arrayToMatrix(gsl_matrix* M , Orbit_Var orb);

// IOmodule.c
void printOrbit(Orbit orb);                   // Prints Orbit struct
void printPoint(Point p);                     // Prints Point struct
void tofileOrbit(FILE * fptr, Orbit orb);
void printArray(int n , double* v);

#endif //_MYHEADER_H
