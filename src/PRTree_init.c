#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rmath.h>  // needed for some random number generators

/* ------------------------------------ */
/* R functions to be passed to FORTRAN  */
/* ------------------------------------ */
/* Distributions */

double F77_SUB(pnormr)(double * x, double * mu, double * sigma){return pnorm(*x, *mu, *sigma, 1,0);} /* Gaussian */

/* --------------------------- */
/* .Fortran calls */
/* --------------------------- */

/* Extracting components */
extern void F77_NAME(pr_treer)(double * y, double * X, int * nrow, int * ncol, double * sigma, int * dim_sigma, int * max_terminal_nodes, double * cp, int * max_depth, int * n_min, double * perc_x, double * p_min, int * Iindep, double * P, int * dim_P, double * gamma, double * yhat, double * MSE, int * nodes_matrix_info, double * cutpoints, double * inf, double * sup, double * sigma_best, int * XRegion);

/* Prediction */
extern void F77_NAME(predict_pr_treer)(int * Iindep, int * nrow, int * ncol, double * X_test, double * inf, double * sup, int * n_terminal_nodes, int * tn, double * P, double * gamma, int * dim_sigma, double * sigma, double * yhat_test);

/* --------------------------- */
/* end of .Fortran calls */
/* --------------------------- */

static const R_FortranMethodDef FortranEntries[] = {
		{"pr_treer",                (DL_FUNC) &F77_NAME(pr_treer),                24},
		{"predict_pr_treer",        (DL_FUNC) &F77_NAME(predict_pr_treer),        13},
    {NULL, NULL, 0}
};

void R_init_PRTree(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
