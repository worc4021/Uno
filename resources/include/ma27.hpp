#ifndef MA27_HPP
#define MA27_HPP
#include "fc_help.h"


extern "C" {
extern void FC_ma27id(  int *ICNTL, double* CNTL);

extern void FC_ma27ad(  int& N, int& NZ, 
                        int* IRN, int* ICN, 
                        int* IW, int& LIW, 
                        int* IKEEP, int* IW1,
                        int& NSTEPS, int& IFLAG,
                        int* ICNTL,double* CNTL, 
                        int* INFO, double& OPS);
    
extern void FC_ma27bd(  int& N, int& NZ, 
                        int* IRN, int* ICN, 
                        double* A, int& LA, 
                        int* IW, int& LIW, 
                        int* IKEEP, int& NSTEPS, 
                        int& MAXFRT, int* IW1, 
                        int* ICNTL, double* CNTL, 
                        int* INFO);
    
extern void  FC_ma27cd( int& N, double* A, 
                        int& LA, int* IW, 
                        int& LIW, double* W, 
                        int& MAXFRT, double* RHS,
                        int* IW1, int& NSTEPS, 
                        int* ICNTL, int* INFO );

}

#endif // MA27_HPP