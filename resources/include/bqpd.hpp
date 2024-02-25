#ifndef BQPD_HPP
#define BQPD_HPP
#include "fc_help.h"

extern "C" {
extern void FC_bqpd(const int* n, const int* m, int* k, int* kmax, double* a, int* la, double* x, double* bl, double* bu, double* f, double* fmin, double* g,
      double* r, double* w, double* e, int* ls, double* alp, int* lp, int* mlp, int* peq, double* ws, int* lws, const int* mode, int* ifail,
      int* info, int* iprint, int* nout);


//  Common
//  ******
//  User information about the lengths of ws and lws is supplied to bqpd in
//    common/wsc/kk,ll,kkk,lll,mxws,mxlws
//  mxws and mxlws must be set to the total length of ws and lws.  kk and ll
//  refer to the length of ws and lws that is used by gdotx. kkk and lll
//  are the numbers of locations used by bqpd and are set by bqpd.

extern struct {
   int kk, ll, kkk, lll, mxws, mxlws;
} FC_wsc;

// fortran common for inertia correction in wdotd
extern struct {
   double alpha;
} FC_kktalphac;

extern struct {
      double eps,tol,emin;
} FC_epsc;

extern struct {
      double sgnf,nrep,npiv,nres;
} FC_repc;

extern struct {
      double nup,nfreq;
} FC_refactorc;

}

#endif // BQPD_HPP