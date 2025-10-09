#include <R.h>
#include <stdlib.h>         // for NULL
#include <R_ext/Rdynload.h> // registering .Fortran
#include <Rmath.h>          // for distribution functions
#include <R_ext/Utils.h>    // for revsort
#include <R_ext/Arith.h>    // for NaN

/* ------------------------------------
  R functions to be passed to FORTRAN
  ------------------------------------ */
/* Sorting */
void F77_SUB(revsortr)(double *x, int *indx, int *n){
    // x*: array to be SORTED (will be MODIFIED)
    // indx*: index array that will be MODIFIED
    // n: size of arrays (will not be modified)
    Rf_revsort(x, indx, *n);
}

/* -------------------------------------------
   Distributions - PURE wrappers for Fortran
   ------------------------------------------ */
/* PURE version of pnorm */
double F77_SUB(pnorm_pure)(double *x, double *mu, double *sigma) {
    // x*: input value pointer
    // mean*: mean parameter pointer
    // sd*: standard deviation pointer
    if (*sigma < 0.0) return R_NaN;
    return Rf_pnorm5(*x, *mu, *sigma, 1, 0);
}

/* PURE version of plnorm */
double F77_SUB(plnorm_pure)(double *x, double *meanlog, double *sdlog) {
    // x*: input value pointer
    // meanlog*: meanlog parameter pointer
    // sdlog*: sdlog parameter pointer
    if (*x <= 0.0) return 0.0;
    if (*sdlog < 0.0) return R_NaN;
    return  Rf_plnorm(*x, *meanlog, *sdlog, 1, 0);
}

/* PURE version of pt */
double F77_SUB(pt_pure)(double *x, double *df) {
    // x*: input value pointer
    // df*: degrees of freedom pointer
    if (*df <= 0.0) return R_NaN;
    return Rf_pt(*x, *df, 1, 0);
}

/* PURE version of pgamma */
double F77_SUB(pgamma_pure)(double *x, double *shape, double *scale) {
    // x*: input value pointer
    // shape*: shape parameter pointer
    // scale*: scale parameter pointer
    if (*x <= 0.0) return 0.0;
    if (*shape < 0.0 || *scale <= 0.0) return R_NaN;
    return Rf_pgamma(*x, *shape, *scale, 1, 0);
}

/* --------------------------- */
/* .Fortran calls */
/* --------------------------- */

/* Build the tree */
extern void F77_NAME(pr_tree_fort)(int *n_obs, int *n_feat, int * n_train, int *idx_train, double *y, double *X, int *n_sigmas, double *sigmas, int *int_param, double *dble_param, int *n_tn, double *P, double *gamma, double *yhat, double *mse, int *nodes, double *thresholds, double *bounds, double *sigma_best, int *Xregion);

/* Prediction */
extern void F77_NAME(predict_pr_tree_fort)(int *dist, double *pdist, int *fill, int *n_obs, int *n_feat, double *X_test, double *bounds, int *n_tn, int *tn, int* nodes_info, double *P, double *gamma, double *sigma, double *yhat_test);

/* --------------------------- */
/* end of .Fortran calls */
/* --------------------------- */

static const R_FortranMethodDef FortranEntries[] = {
    {"pr_tree_fort", (DL_FUNC)&F77_NAME(pr_tree_fort), 20},
    {"predict_pr_tree_fort", (DL_FUNC)&F77_NAME(predict_pr_tree_fort), 14},
    {NULL, NULL, 0}};

void R_init_PRTree(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
