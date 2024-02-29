#ifndef MA27_HPP
#define MA27_HPP
#include "fc_help.h"


extern "C" {
extern void FC_ma27id(int icntl[30], double cntl[5]);

extern void FC_ma27ad(int* n, int* nz, int* irn, int* icn, int iw[],
                        int* liw, int* ikeep, int* iw1, int* nsteps,
                        int* iflag, int icntl[30], double cntl[5],
                        int info[20], double* ops );
    
extern void FC_ma27bd(int* n, int* nz, int* irn, int* icn,
                   double* a, int* la, int iw[], int* liw,
                   int* ikeep, int* nsteps, int* maxfrt, int* iw1,
                   int icntl[30], double cntl[5], int info[20]);
    
extern void  FC_ma27cd(int* n, double* a, int* la, int iw[], int* liw,
                       double* w, int* maxfrt, double* rhs, int* iw1,
                       int* nsteps, int icntl[30], int info[20]);

extern void FC_ma27dd(  int* job, int* n, 
                        int* ne,  double* a, 
                        int* irn, int* jcn, 
                        double* fact, int* lfact,
                        int* ifact, int* lifact, 
                        double* rhs, double* x, 
                        double* resid, double* work,
                        int* iwork, int* icntl,
                        double* cntl, int* info,
                        double* rinfo);
}

#endif // MA27_HPP