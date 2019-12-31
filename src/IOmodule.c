

#include "../myHeader.h"

// Print in terminal
void printOrbit(Orbit orb){
    for (int i=0;i<orb.nt;++i){
        printf("%g %g %g %g %g\n",orb.t[i],orb.x1[i],orb.x2[i],orb.x3[i],orb.x4[i]);
    }
}

void printPoint(Point p){
    printf("%g %g %g %g %g\n",p.t,p.x[0],p.x[1],p.x[2],p.x[3]);
}

void printArray(int n , double* v){
    for (int i = 0 ; i < n ; ++i) printf("%g, " , v[i]);
    printf("\n");
}


// Print to file
void tofileOrbit(FILE * fptr, Orbit orb){
    fprintf(fptr , "%d\n", orb.nt);
    for (int i=0;i<orb.nt;++i){
        fprintf(fptr , "%12.6g %12.6g %12.6g %12.6g %12.6g \n",orb.t[i],orb.x1[i],orb.x2[i],orb.x3[i],orb.x4[i]);
    }
}
