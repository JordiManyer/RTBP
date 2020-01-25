#include "../myHeader.h"

// -----------------------------------------------------------------------------------------------------
// -----------------------------------------------------------------------------------------------------
/* Computes the Poincaré section up to 'nCrossings' of 'func',
 * starting from ('tstart', 'x0') in direction 'dir'.  */
Point computePoincare(double* x0 , double tstart , int dir , double (*event)(double*), int nCrossings){
    int itmp, order;
    double atol, rtol;
    MY_FLOAT  startT, stopT, nextT;
    MY_FLOAT  xx[5];

    int numC;
    double tlast, tnew, tstep;
    double flast,fnew;
    double xlast[4],xnew[4];
    Point res;

    // Tolerances
    itmp = 0;
    order = 20;
    atol = log10(1.e-18);
    rtol = log10(1.e-16);

    // Initial conditions
    numC = 0;
    if (dir == 1) {tstep = 0.1;} else {tstep = -0.1;}
    for (int i = 0 ; i < 4 ; ++i) xlast[i] = x0[i];
    tlast = tstart;
    flast = event(xlast);

    // First iterate all the first crossings fast
    while (numC < nCrossings-1 && tlast < 1.e4){
        tnew = tlast + tstep; // integrate up to here

        // Integration times
        MakeMyFloatA(startT, tlast);
        MakeMyFloatA(stopT , tnew);
        MakeMyFloatA(nextT , tstep);

        for (int i = 0 ; i < 4 ; ++i) xx[i] = xlast[i];
        do  { // Integrate up to tnew
            itmp = taylor_step_rtbp( &startT, xx, dir, 2, atol, rtol, &stopT, &nextT, &order);
            if (verbose > 1) printf("%g %g %g %g %g\n",(double)xx[0],(double)xx[1],(double)xx[2],(double)xx[3],(double)startT);
        } while(itmp == 0);
        for (int i = 0 ; i < 4 ; ++i) xnew[i] = xx[i]; tnew = startT;

        // If we cross the poincare section, increase numC
        fnew = event(xnew);
        if(flast*fnew < 0) numC++;

        // Update everything else regardless
        tlast = tnew;
        flast = fnew;
        for (int i = 0 ; i < 4 ; ++i) xlast[i] = xnew[i];
    }

    // The last crossing needs to be treated carefully:
    if (dir == 1) {tstep = 0.1;} else {tstep = -0.1;}
    while (fabs(tstep) > PCTol && tlast < 1.e4){ // While not having a certain tolerance
        tnew = tlast + tstep; // integrate up to here

        // Integration times
        MakeMyFloatA(startT, tlast);
        MakeMyFloatA(stopT , tnew);
        MakeMyFloatA(nextT , tstep);

        for (int i = 0 ; i < 4 ; ++i) xx[i] = xlast[i];
        do  { // Integrate up to tnew
            itmp = taylor_step_rtbp( &startT, xx, dir, 2, atol, rtol, &stopT, &nextT, &order);
            if (verbose > 1) printf("%g \n",tstep);
            if (verbose > 1) printf("%g %g %g %g %g\n",(double)xx[0],(double)xx[1],(double)xx[2],(double)xx[3],(double)startT);
        } while(itmp == 0);
        for (int i = 0 ; i < 4 ; ++i) xnew[i] = xx[i]; tnew = startT;

        fnew = event(xnew);
        if(flast*fnew < 0) {
            // If we cross the poincare section, we decrease tstep but do not advance.
            tstep = tstep / 10.0;
        } else {
            // Otherwise we advance
            tlast = tnew;
            flast = fnew;
            for (int i = 0 ; i < 4 ; ++i) xlast[i] = xnew[i];
        }
    }

    res.t = tnew;
    for (int i = 0 ; i < 4 ; ++i) res.x[i] = xnew[i];
    return res;
}


void computeOrbitPoincare(Orbit orb, double* x0 , double tstart , int nPoints , int dir , double (*event)(double*)) {
    Point p;
    double x[4]; double time;

    time = tstart; x[0] = x0[0]; x[1] = x0[1]; x[2] = x0[2]; x[3] = x0[3]; // Initialize buffers
    orb.t[0] = tstart; orb.x1[0] = x[0]; orb.x2[0] = x[1]; orb.x3[0] = x[2]; orb.x4[0] = x[3]; // Save first point
    for (int i = 1 ; i < nPoints ; ++i) {
        p = computePoincare(x , time , dir , event, 1); // Compute Poincare map step

        // Save values into orbit
        orb.t[i] = p.t; orb.x1[i] = p.x[0]; orb.x2[i] = p.x[1]; orb.x3[i] = p.x[2]; orb.x4[i] = p.x[3];
        // Update buffers
        time = p.t; x[0] = p.x[0]; x[1] = p.x[1]; x[2] = p.x[2]; x[3] = p.x[3];
    }
    // Returns the orbit
}

