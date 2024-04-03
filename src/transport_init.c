#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
// consider including RcppExports.h here (after renaming)
// somehow the current solution with two files is a bit awkward
// and hard to maintain


/* .C calls */
// 221017: I think function package_native_routine_registration_skeleton is still producing a legit result
// with all the void * arguments, except that it produces e.g. aha_free(), which in C is a
// function with unspecified arguments (rather than a function with no arguments) and therefore
// leads to the function without prototype warning (not allowed in C2x)
extern void aha_compute_transport(int *n, double *x, double *y, double *w, double *source_measure, int *res);
extern void aha_dphi(int *n, double *x, double *y, double *w, double *source_measure,
                     double *target_measure, int *exact, double *res);
extern void aha_free(void);  // or void * in this case as functions here should only take pointers
                             // on the other hand void is not a variable either but just a "keyword", no?
extern void aha_get_transport(int *size, double *from, double *to, double *mass);
extern void aha_init(int *n, int *m, double *rect, int *npoints);
extern void aha_phi(int *n, double *x, double *y, double *w, double *source_measure, double *target_measure,
                    int *exact, double *res);
extern void aha_wasserstein(int *n, double *x, double *y, double *w, double *source_measure, double *res);
extern void auction(int *desirem, int *nn, int *pers_to_obj, double *price, int *kk, double *eps);
extern void auctionbf(int *desirem, int *nn, int *pers_to_obj, double *price, double *profit, int *kk, double *eps);
extern void compute_power_diagram(int *cell_size, int *n, double *x, double *y, double *w, double *rect);
extern void decompose_c(int *n, double *x, double *y, double *m, int *n0, double *x0, double *y0,
                             double *m0, int *p, double *eps);
extern void get_power_diagram(int *size, double *x, double *y);
extern void primaldual(int *d, int *rmass, int *cmass, int *numr, int *numc, int *flowmatrix);
extern void revsimplex(int *mm, int *nn, int *a, int *b, double *costm,
                       int *assignment, int *basis, int *startgiven);
extern void shortsimplex(int *ss, int *kk, double *pp, int *mm, int *nn, int *a, int *b,
                         double *costm, int *assignment, int *basis);

/* .Call calls */
// 221017: this I copied originally from RcppExports.cpp I pressume (surprisingly Rcpp::compileAttributes
// seems to recognize this and does not generate it again). Since C-functions without prototype
// are not allowed anymore, and fun() is not regarded as a prototype (means unspecified arguments
// in C but no arguments in C++), I added void here. According to forums fun(void) should be
// exactly the same as fun() in C++, so the literal discrepancy with the actual definitions 
// should not matter (some are in RcppExports, so I cannot / do not want to change it there).
extern SEXP _transport_create_diagram(SEXP);
extern SEXP _transport_cgal_present(void);
extern SEXP _transport_semidiscrete_p1(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _transport_cplex_present(void);
extern SEXP _transport_SolveHierarchicalTransport(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _transport_networkflow(SEXP, SEXP, SEXP, SEXP);
extern SEXP _transport_gen_cost(SEXP, SEXP, SEXP);
extern SEXP _transport_openmp_present(void);

static const R_CMethodDef CEntries[] = {
    {"aha_compute_transport", (DL_FUNC) &aha_compute_transport,  6},
    {"aha_dphi",              (DL_FUNC) &aha_dphi,               8},
    {"aha_free",              (DL_FUNC) &aha_free,               0},
    {"aha_get_transport",     (DL_FUNC) &aha_get_transport,      4},
    {"aha_init",              (DL_FUNC) &aha_init,               4},
    {"aha_phi",               (DL_FUNC) &aha_phi,                8},
    {"aha_wasserstein",       (DL_FUNC) &aha_wasserstein,        6},
    {"auction",               (DL_FUNC) &auction,                6},
    {"auctionbf",             (DL_FUNC) &auctionbf,              7},
    {"compute_power_diagram", (DL_FUNC) &compute_power_diagram,  6},
    {"decompose_c",           (DL_FUNC) &decompose_c,           10},
    {"get_power_diagram",     (DL_FUNC) &get_power_diagram,      3},
    {"primaldual",            (DL_FUNC) &primaldual,             6},
    {"revsimplex",            (DL_FUNC) &revsimplex,             8},
    {"shortsimplex",          (DL_FUNC) &shortsimplex,          10},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"_transport_create_diagram", (DL_FUNC) &_transport_create_diagram, 1},
    {"_transport_cgal_present", (DL_FUNC) &_transport_cgal_present,0},
    {"_transport_semidiscrete_p1", (DL_FUNC) &_transport_semidiscrete_p1, 6},
    {"_transport_cplex_present", (DL_FUNC) &_transport_cplex_present, 0},
    {"_transport_SolveHierarchicalTransport", (DL_FUNC) &_transport_SolveHierarchicalTransport, 13},
    {"_transport_networkflow", (DL_FUNC) &_transport_networkflow, 4},
    {"_transport_gen_cost", (DL_FUNC) &_transport_gen_cost, 3},
    {"_transport_openmp_present", (DL_FUNC) &_transport_openmp_present, 0},
    {NULL, NULL, 0}
};

void R_init_transport(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
