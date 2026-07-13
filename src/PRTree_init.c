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
   Distributions - Vectorized wrappers for Fortran
   ------------------------------------------ */
/* Vectorized version of pnorm */
void F77_SUB(pnorm_v)(int *n, double *v, double *mu, double *sigma, double *out) {
    for (int i = 0; i < *n; i++) {
        out[i] = Rf_pnorm5(*v, mu[i], *sigma, 1, 0);
    }
}

/* Vectorized version of plnorm */
void F77_SUB(plnorm_v)(int *n, double *v, double *mu, double *sigma, double *sdlog, double *out) {
    double inv_sigma = 1.0 / (*sigma);
    for (int i = 0; i < *n; i++) {
        double z = (*v - mu[i]) * inv_sigma;
        out[i] = Rf_plnorm(z, 0.0, *sdlog, 1, 0);
    }
}

/* Vectorized version of pt */
void F77_SUB(pt_v)(int *n, double *v, double *mu, double *sigma, double *df, double *out) {
    double inv_sigma = 1.0 / (*sigma);
    for (int i = 0; i < *n; i++) {
        double z = (*v - mu[i]) * inv_sigma;
        out[i] = Rf_pt(z, *df, 1, 0);
    }
}

/* Vectorized version of pgamma */
void F77_SUB(pgamma_v)(int *n, double *v, double *mu, double *sigma, double *shape, double *out) {
    double inv_sigma = 1.0 / (*sigma);
    for (int i = 0; i < *n; i++) {
        double z = (*v - mu[i]) * inv_sigma;
        out[i] = Rf_pgamma(z, *shape, 1.0, 1, 0);
    }
}

/* --------------------------- */
/* .Fortran calls */
/* --------------------------- */

/* Build the tree */
extern void F77_NAME(pr_tree_fort)(int *n_obs, int *n_feat, int *n_train, double *y_train, double *y_test, double *x_train, double *x_test, int *n_sigmas, double *sigmas, int *int_param, double *dble_param, int *n_tn, double *p_train, double *p_test, double *gammahat, double *yhat_train, double *yhat_test, double *mse, int *nodes_info, double *thresholds, double *sigma_best, int *xregion);

/* Prediction */
extern void F77_NAME(predict_pr_tree_fort)(int *dist, double *pdist, int *fill, int *n_obs, int *n_feat, double *X_test, double *thresholds, int *n_tn, int *tn, int* nodes_info, double *P, double *gamma, double *sigma, double *yhat_test);

/* --------------------------- */
/* end of .Fortran calls */
/* --------------------------- */

static const R_FortranMethodDef FortranEntries[] = {
    {"pr_tree_fort", (DL_FUNC)&F77_NAME(pr_tree_fort), 22},
    {"predict_pr_tree_fort", (DL_FUNC)&F77_NAME(predict_pr_tree_fort), 14},
    {NULL, NULL, 0}};

void R_init_PRTree(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
