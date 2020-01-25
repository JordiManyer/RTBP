#include "../myHeader.h"

/* Computes the orbit starting from ('tstart', 'x0'),
 * in direction 'dir' and ending at 'tend'.  */
Orbit computeOrbit(double* x0 , double tstart , double tend , double tstep , int dir){
    int itmp, order;
    double atol, rtol;
    MY_FLOAT  startT, stopT, nextT;
    MY_FLOAT  xx[5];

    // Orbit declaration
    Orbit res;
    double t[orbitMaxLength], x1[orbitMaxLength], x2[orbitMaxLength], x3[orbitMaxLength], x4[orbitMaxLength];
    res.t = t; res.x1 = x1; res.x2 = x2; res.x3 = x3; res.x4 = x4;

    itmp = 0;
    order = 20;

    // Initial conditions
    MakeMyFloatA(xx[0] ,x0[0]);
    MakeMyFloatA(xx[1] ,x0[1]);
    MakeMyFloatA(xx[2] ,x0[2]);
    MakeMyFloatA(xx[3] ,x0[3]);

    // Integration times
    MakeMyFloatA(startT, tstart);
    MakeMyFloatA(stopT , tend);
    MakeMyFloatA(nextT , tstep);

    // Tolerances
    atol = log10(1.e-18);
    rtol = log10(1.e-16);

    res.nt = 1;
    res.t[0]  = tstart;
    res.x1[0] = x0[0];
    res.x2[0] = x0[1];
    res.x3[0] = x0[2];
    res.x4[0] = x0[3];
    if(tstep < (double)0.0) { dir = -1;}
    do  {
        itmp = taylor_step_rtbp( &startT, xx, dir, 2, atol, rtol, &stopT, &nextT, &order);
        res.t[res.nt] = startT; res.x1[res.nt] = xx[0];
        res.x2[res.nt] = xx[1]; res.x3[res.nt] = xx[2];
        res.x4[res.nt] = xx[3];
        res.nt++;
        if (verbose > 1) printf("%g %g %g %g %g\n",(double)xx[0],(double)xx[1],(double)xx[2],(double)xx[3],(double)startT);
    } while(itmp == 0);

    return res;
}

// -----------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------
void computeOrbitNpoints(Orbit orb, double* x0 , double tstart , double tend , int nPoints , int dir){
    int itmp, order;
    double atol, rtol;
    MY_FLOAT  startT, stopT, nextT;
    MY_FLOAT  xx[5];

    double tlast, tnew, tstep;
    double xlast[4],xnew[4];

    // Tolerances
    itmp = 0;
    order = 20;
    atol = log10(1.e-18);
    rtol = log10(1.e-16);

    // Initial conditions
    tstep = (tend-tstart)/(nPoints-1.0);
    for (int i = 0 ; i < 4 ; ++i) xlast[i] = x0[i];
    orb.t[0] = tstart; orb.x1[0] = x0[0]; orb.x2[0] = x0[1]; orb.x3[0] = x0[2]; orb.x4[0] = x0[3];
    tlast = tstart;

    // First iterate all the first crossings fast
    for (int j = 1; j < nPoints ; ++j){
        tnew = tlast + tstep; // integrate up to here

        // Integration times
        MakeMyFloatA(startT, tlast);
        MakeMyFloatA(stopT , tnew);
        MakeMyFloatA(nextT , tstep);

        for (int i = 0 ; i < 4 ; ++i) xx[i] = xlast[i];
        do  { // Integrate up to tnew
            itmp = taylor_step_rtbp( &startT, xx, dir, 1, atol, rtol, &stopT, &nextT, &order);
            if (verbose > 1) printf("%g %g %g %g %g\n",(double)xx[0],(double)xx[1],(double)xx[2],(double)xx[3],(double)startT);
        } while(itmp == 0);
        for (int i = 0 ; i < 4 ; ++i) xnew[i] = xx[i]; tnew = startT;
        // Update everything
        tlast = tnew;
        for (int i = 0 ; i < 4 ; ++i) xlast[i] = xnew[i];
        orb.t[j] = tnew; orb.x1[j] = xnew[0]; orb.x2[j] = xnew[1]; orb.x3[j] = xnew[2]; orb.x4[j] = xnew[3];
    }
}



// -----------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------
void computeVariationalsNpoints(Orbit_Var orb, double* x0 , double tstart , double tend , int nPoints , int dir){
    int itmp, order;
    double atol, rtol;
    MY_FLOAT  startT, stopT, nextT;
    MY_FLOAT  xx[21];

    double tlast, tnew, tstep;
    double xlast[20],xnew[20];

    // Tolerances
    itmp = 0;
    order = 20;
    atol = log10(1.e-18);
    rtol = log10(1.e-16);

    // Initial conditions
    tstep = (tend-tstart)/(nPoints-1.0);
    for (int i = 0 ; i < 20 ; ++i) xlast[i] = x0[i];
    orb.t[0] = tstart;
    orb.x1[0] = x0[0]; orb.x2[0] = x0[1]; orb.x3[0] = x0[2]; orb.x4[0] = x0[3];

    orb.x5[0] = x0[4];   orb.x6[0] = x0[5];   orb.x7[0] = x0[6];   orb.x8[0] = x0[7];
    orb.x9[0] = x0[8];   orb.x10[0] = x0[9];  orb.x11[0] = x0[10]; orb.x12[0] = x0[11];
    orb.x13[0] = x0[12]; orb.x14[0] = x0[13]; orb.x15[0] = x0[14]; orb.x16[0] = x0[15];
    orb.x17[0] = x0[16]; orb.x18[0] = x0[17]; orb.x19[0] = x0[18]; orb.x20[0] = x0[19];
    tlast = tstart;

    // First iterate all the first crossings fast
    for (int j = 1; j < nPoints ; ++j){
        tnew = tlast + tstep; // integrate up to here

        // Integration times
        MakeMyFloatA(startT, tlast);
        MakeMyFloatA(stopT , tnew);
        MakeMyFloatA(nextT , tstep);

        for (int i = 0 ; i < 20 ; ++i) xx[i] = xlast[i];
        do  { // Integrate up to tnew
            itmp = taylor_step_rtbp_variationals( &startT, xx, dir, 1, atol, rtol, &stopT, &nextT, &order);
            if (verbose > 1) printf("%g %g %g %g %g\n",(double)xx[0],(double)xx[1],(double)xx[2],(double)xx[3],(double)startT);
        } while(itmp == 0);
        for (int i = 0 ; i < 20 ; ++i) xnew[i] = xx[i]; tnew = startT;
        // Update everything
        tlast = tnew;
        for (int i = 0 ; i < 20 ; ++i) xlast[i] = xnew[i];
        orb.t[j] = tnew; orb.x1[j] = xnew[0]; orb.x2[j] = xnew[1]; orb.x3[j] = xnew[2]; orb.x4[j] = xnew[3];

        orb.x5[j]  = xnew[4];  orb.x6[j]  = xnew[5];   orb.x7[j] = xnew[6];   orb.x8[j] = xnew[7];
        orb.x9[j]  = xnew[8];  orb.x10[j] = xnew[9];  orb.x11[j] = xnew[10]; orb.x12[j] = xnew[11];
        orb.x13[j] = xnew[12]; orb.x14[j] = xnew[13]; orb.x15[j] = xnew[14]; orb.x16[j] = xnew[15];
        orb.x17[j] = xnew[16]; orb.x18[j] = xnew[17]; orb.x19[j] = xnew[18]; orb.x20[j] = xnew[19];
    }
}