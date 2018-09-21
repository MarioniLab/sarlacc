#include "umi_clusterer.h"

Rcpp::List umi_clusterer::cluster() {
    const size_t n_stored=storage.size();
    if (n_stored > ordering.size()) { 
        ordering.resize(n_stored);
        remaining.resize(n_stored);
    }
    std::iota(ordering.begin(), ordering.begin() + n_stored, 0);

    // Removing solo reads beforehand.
    size_t infront=0, n_out=0;
    for (size_t a=0; a<n_stored; ++a) {
        const auto& curlen=(remaining[a]=storage.get_len(a));
        if (curlen > 1) {
            continue;
        } else if (curlen==1) {
            if (static_cast<size_t>(*storage.get_start(a))!=a) {
                throw std::runtime_error("single-read groups should contain only the read itself");
            }

            Rcpp::IntegerVector tosave(1, a);
            if (n_out < output.size()) {
                output[n_out]=tosave;
            } else {
                output.push_back(tosave);
            }
            ++n_out;
        }

        std::swap(ordering[a], ordering[infront]);
        ++infront;
    }

    // Keeping the top nodes, accounting for updated counts after processing prior nodes.
    auto left=ordering.begin() + infront, right=ordering.begin() + n_stored;
    while (left < right) {

        // Wiping out empty nodes.
        auto discarded=left;
        while (left!=right && remaining[*left]==0) {
            std::swap(*left, *discarded);
            ++discarded;
            ++left;
        }
        if (left==right) { 
            break;
        }

        // Finding the node with the maximum size.
        auto maxIt=std::max_element(left, right,
            [&] (const size_t& L, const size_t& R) -> bool {
                if (remaining[L]==remaining[R]) {
                    return L < R;
                }
                return remaining[L] < remaining[R];
            }
        );

        const auto curstart=storage.get_start(*maxIt);
        const auto curend=curstart + storage.get_len(*maxIt);
        --right;
        std::swap(*maxIt, *right); // swapping it out of [left, right), effectively a pop_back().

        // Adding neighbors if they have not already been used somewhere else.
        size_t n_cluster=0;
        for (auto curIt=curstart; curIt<curend; ++curIt) { 
            const auto& neighbor=*curIt;
            auto& remain=remaining[neighbor];
            if (remain==0) {
                continue;
            }

            if (n_cluster < per_cluster_values.size()) {
                per_cluster_values[n_cluster]=neighbor;
            } else {
                per_cluster_values.push_back(neighbor);
            }
            ++n_cluster;
            remain=0;

            // Also decrementing the counts of all neighbors of neighbors.
            const auto neighstart=storage.get_start(neighbor);
            const auto neighend=neighstart + storage.get_len(neighbor);
            for (auto nextIt=neighstart; nextIt<neighend; ++nextIt) {
                const auto& nextneighbor=*nextIt;
                auto& nextremain=remaining[nextneighbor];
                if (nextremain > 0) { 
                    --nextremain;
                }
            }
        }

        if (n_cluster) {
            Rcpp::IntegerVector tosave(per_cluster_values.begin(), per_cluster_values.begin() + n_cluster);
            if (n_out < output.size()) {
                output[n_out]=tosave;
            } else {
                output.push_back(tosave);
            }
            ++n_out;
        }
    }

    return Rcpp::List(output.begin(), output.end());
}

/*****************************
 * R-level testing function. *
 *****************************/

SEXP cluster_umis_test (SEXP links) {
    BEGIN_RCPP
    Rcpp::List Links(links);
    const size_t nsets=Links.size();
    umi_clusterer clust;

    for (size_t l=0; l<nsets; ++l) {
        Rcpp::IntegerVector current=Links[l];
        clust.storage.add(current.begin(), current.end());
        auto startIt=clust.storage.get_start_unsafe(l);
        for (size_t i=0; i<current.size(); ++i, ++startIt) {
            --(*startIt); // get to zero indexing.
        }
    }

    Rcpp::List curout=clust.cluster();
    for (size_t i=0; i<curout.size(); ++i) {
        Rcpp::IntegerVector curvec=curout[i];
        for (auto& x : curvec) { --x; }
        curout[i]=curvec;
    }
    return curout;
    END_RCPP
}
