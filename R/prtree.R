##' @title
##' PRTree: Probabilistic Regression Tress
##'
##' @description
##' Probabilistic Regression Trees (PRTree). Functions for fitting and predicting
##' PRTree models with some adaptations to handle missing values. The main
##' calculations are performed in 'FORTRAN', resulting in highly efficient
##' algorithms. This package's implementation is based on the PRTree methodology
##' described in Alkhoury, S.; Devijver, E.; Clausel, M.; Tami, M.; Gaussier, E.;
##' Oppenheim, G. (2020) - "Smooth And Consistent Probabilistic Regression Trees"
##' <https://proceedings.neurips.cc/paper_files/paper/2020/file/8289889263db4a40463e3f358bb7c7a1-Paper.pdf>.
##'
##' @author Taiane Schaedler Prass \email{taianeprass@@gmail.com} and Alisson Silva Neimaier \email{alissonneimaier@hotmail.com}
##'
##' @docType package
##'
##' @name PRTree.Package
##'
##' @aliases PRTree
##'
##' @keywords internal
"_PACKAGE"
##' @useDynLib PRTree, .registration=TRUE
##'
NULL
## > NULL
