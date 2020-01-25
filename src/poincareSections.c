
#include "../myHeader.h"

double xequalzero(double *x){
    return x[0];
}

double yequalzero(double *x){
    return x[1];
}

double xequalvalue(double *x){
    return x[0] - PoincareValue;
}

double xprimeequalzero(double *x){
    return x[2];
}
