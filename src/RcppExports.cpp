// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// create_diagram
NumericMatrix create_diagram(NumericMatrix sites);
RcppExport SEXP _transport_create_diagram(SEXP sitesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type sites(sitesSEXP);
    rcpp_result_gen = Rcpp::wrap(create_diagram(sites));
    return rcpp_result_gen;
END_RCPP
}
// gen_cost
NumericMatrix gen_cost(NumericMatrix AR, NumericMatrix BR, int threads);
RcppExport SEXP _transport_gen_cost(SEXP ARSEXP, SEXP BRSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type AR(ARSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type BR(BRSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(gen_cost(AR, BR, threads));
    return rcpp_result_gen;
END_RCPP
}
// gen_cost0
NumericMatrix gen_cost0(NumericMatrix xx, NumericMatrix yy);
RcppExport SEXP _transport_gen_cost0(SEXP xxSEXP, SEXP yySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type yy(yySEXP);
    rcpp_result_gen = Rcpp::wrap(gen_cost0(xx, yy));
    return rcpp_result_gen;
END_RCPP
}
// gen_cost0d
NumericMatrix gen_cost0d(NumericMatrix xx, NumericMatrix yy);
RcppExport SEXP _transport_gen_cost0d(SEXP xxSEXP, SEXP yySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type yy(yySEXP);
    rcpp_result_gen = Rcpp::wrap(gen_cost0d(xx, yy));
    return rcpp_result_gen;
END_RCPP
}
// networkflow
List networkflow(NumericMatrix a, NumericMatrix b, NumericMatrix C, int threads, size_t maxiters);
RcppExport SEXP _transport_networkflow(SEXP aSEXP, SEXP bSEXP, SEXP CSEXP, SEXP threadsSEXP, SEXP maxitersSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type a(aSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type b(bSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type C(CSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< size_t >::type maxiters(maxitersSEXP);
    rcpp_result_gen = Rcpp::wrap(networkflow(a, b, C, threads, maxiters));
    return rcpp_result_gen;
END_RCPP
}
// openmp_present
SEXP openmp_present();
RcppExport SEXP _transport_openmp_present() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(openmp_present());
    return rcpp_result_gen;
END_RCPP
}
// cplex_present
SEXP cplex_present();
RcppExport SEXP _transport_cplex_present() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(cplex_present());
    return rcpp_result_gen;
END_RCPP
}
// cgal_present
SEXP cgal_present();
RcppExport SEXP _transport_cgal_present() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(cgal_present());
    return rcpp_result_gen;
END_RCPP
}
// semidiscrete_p1
List semidiscrete_p1(NumericMatrix source_mat, NumericMatrix target_mat, bool verbose, bool target_in_genpos, double regularization_strength, NumericMatrix transportplan);
RcppExport SEXP _transport_semidiscrete_p1(SEXP source_matSEXP, SEXP target_matSEXP, SEXP verboseSEXP, SEXP target_in_genposSEXP, SEXP regularization_strengthSEXP, SEXP transportplanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type source_mat(source_matSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type target_mat(target_matSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type target_in_genpos(target_in_genposSEXP);
    Rcpp::traits::input_parameter< double >::type regularization_strength(regularization_strengthSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type transportplan(transportplanSEXP);
    rcpp_result_gen = Rcpp::wrap(semidiscrete_p1(source_mat, target_mat, verbose, target_in_genpos, regularization_strength, transportplan));
    return rcpp_result_gen;
END_RCPP
}
// SolveHierarchicalTransport
SEXP SolveHierarchicalTransport(SEXP x, SEXP y, SEXP xydepth, SEXP xydimensions, SEXP compdepth, SEXP measureScaleVecPre, SEXP keepBasisVecPre, SEXP refineBasisVecPre, SEXP layerCoarsestVecPre, SEXP verboseVecPre, SEXP assignment, SEXP udummy, SEXP vdummy);
RcppExport SEXP _transport_SolveHierarchicalTransport(SEXP xSEXP, SEXP ySEXP, SEXP xydepthSEXP, SEXP xydimensionsSEXP, SEXP compdepthSEXP, SEXP measureScaleVecPreSEXP, SEXP keepBasisVecPreSEXP, SEXP refineBasisVecPreSEXP, SEXP layerCoarsestVecPreSEXP, SEXP verboseVecPreSEXP, SEXP assignmentSEXP, SEXP udummySEXP, SEXP vdummySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type y(ySEXP);
    Rcpp::traits::input_parameter< SEXP >::type xydepth(xydepthSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xydimensions(xydimensionsSEXP);
    Rcpp::traits::input_parameter< SEXP >::type compdepth(compdepthSEXP);
    Rcpp::traits::input_parameter< SEXP >::type measureScaleVecPre(measureScaleVecPreSEXP);
    Rcpp::traits::input_parameter< SEXP >::type keepBasisVecPre(keepBasisVecPreSEXP);
    Rcpp::traits::input_parameter< SEXP >::type refineBasisVecPre(refineBasisVecPreSEXP);
    Rcpp::traits::input_parameter< SEXP >::type layerCoarsestVecPre(layerCoarsestVecPreSEXP);
    Rcpp::traits::input_parameter< SEXP >::type verboseVecPre(verboseVecPreSEXP);
    Rcpp::traits::input_parameter< SEXP >::type assignment(assignmentSEXP);
    Rcpp::traits::input_parameter< SEXP >::type udummy(udummySEXP);
    Rcpp::traits::input_parameter< SEXP >::type vdummy(vdummySEXP);
    rcpp_result_gen = Rcpp::wrap(SolveHierarchicalTransport(x, y, xydepth, xydimensions, compdepth, measureScaleVecPre, keepBasisVecPre, refineBasisVecPre, layerCoarsestVecPre, verboseVecPre, assignment, udummy, vdummy));
    return rcpp_result_gen;
END_RCPP
}
