#ifndef CLUSTER_UMIS_H
#define CLUSTER_UMIS_H

#include "Rcpp.h"
#include "value_store.h"

Rcpp::List cluster_umis(const value_store<int>&);

#endif
