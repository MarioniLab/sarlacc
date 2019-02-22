#include "sarlacc.h"
#include "cluster_umis.h"

/*****************************
 * R-level testing function. *
 *****************************/

SEXP cluster_umis_test (SEXP links) {
    BEGIN_RCPP
    Rcpp::List Links(links);
    const size_t nsets=Links.size();
    value_store<int> storage;

    for (size_t l=0; l<nsets; ++l) {
        Rcpp::IntegerVector current=Links[l];
        storage.add(current.begin(), current.end());
        auto startIt=storage.get_start(l);
        for (size_t i=0; i<current.size(); ++i, ++startIt) {
            --(*startIt); // get to zero indexing.
        }
    }

    Rcpp::List curout=cluster_umis(storage);
    for (size_t i=0; i<curout.size(); ++i) {
        Rcpp::IntegerVector curvec=curout[i];
        for (auto& x : curvec) { ++x; } // get back to 1-indexing.
        curout[i]=curvec;
    }
    return curout;
    END_RCPP
}
