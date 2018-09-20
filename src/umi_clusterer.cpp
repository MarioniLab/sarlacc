#include "umi_clusterer.h"

umi_clusterer::umi_clusterer () : n_stored(0), nvals_stored(0), storage_values(1000), storage_start(100), storage_len(100) {}

void umi_clusterer::clear() { 
    n_stored=0; 
    nvals_stored=0; 
    return;
}

Rcpp::List umi_clusterer::cluster() {
    if (n_stored > ordering.size()) { 
        ordering.resize(n_stored);
        remaining.resize(n_stored);
    }
    std::iota(ordering.begin(), ordering.begin() + n_stored, 0);
    std::copy(storage_len.begin(), storage_len.begin() + n_stored, remaining.begin());

    // Removing solo reads beforehand.
    size_t infront=0;
    for (size_t a=0; a<n_stored; ++a) {
        const auto& curlen=storage_len[a];
        if (curlen!=1) {
            continue;
        }
        
        if (static_cast<size_t>(storage_values[storage_start[a]])!=a) {
            throw std::runtime_error("single-read groups should contain only the read itself");
        }
        std::swap(ordering[a], ordering[infront]);
        ++infront;
    }

    // Keeping the top nodes, accounting for updated counts after processing prior nodes.
    size_t left=infront, right=n_stored, n_out=0;
    while (left < right) {

        // Wiping out empty nodes.
        size_t discarded=left;
        while (left < right) {
            const size_t& current=remaining[left];
            if (current==0) {
                std::swap(ordering[left], ordering[discarded]);
                ++discarded;
            } 
            ++left;
        }
        if (left==right) { 
            break;
        }

        // Finding the node with the maximum size.
        auto maxIt=std::max_element(ordering.begin(), ordering.end(), 
            [&] (const size_t& left, const size_t& right) -> bool {
                if (remaining[left]==remaining[right]) {
                    return left < right;
                }
                return remaining[left] < remaining[right];
            }
        );

        const auto curstart=storage_start[*maxIt];
        const auto curend=curstart + storage_len[*maxIt];
        std::swap(*maxIt, ordering[right-1]);
        --right;

        // Adding neighbors if they have not already been used somewhere else.
        size_t n_cluster=0;        
        for (auto i=curstart; i<curend; ++i) {
            const auto neighbor=storage_values[i];
            auto& remain=remaining[neighbor];
            if (remain==0) {
                continue;
            }

            if (n_cluster >= per_cluster_values.size()) {
                per_cluster_values.push_back(neighbor);
            } else {
                per_cluster_values[n_cluster]=neighbor;
            }
            ++n_cluster;
            remain=0;

            // Also decrementing the counts of all neighbors of neighbors.
            const auto neighstart=storage_start[neighbor];
            const auto neighend=neighstart + storage_len[neighbor];
            for (auto n=neighstart; n<neighend; ++n) {
                const auto nextneighbor=storage_values[n];
                auto& nextremain=remaining[nextneighbor];
                if (nextremain > 0) { 
                    --nextremain;
                }
            }
        }
                
        Rcpp::IntegerVector tosave(per_cluster_values.begin(), per_cluster_values.begin() + n_cluster);
        if (n_out >= output.size()) {
            output[n_out]=tosave;
        } else {
            output.push_back(tosave);
        }
        ++n_out;
    }
    
    return Rcpp::List(output.begin(), output.end());
}
