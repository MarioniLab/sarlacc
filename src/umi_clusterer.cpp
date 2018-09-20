#include "umi_clusterer.h"

Rcpp::List umi_clusterer::cluster() {
    const size_t n_stored=storage.size();
    if (n_stored > ordering.size()) { 
        ordering.resize(n_stored);
        remaining.resize(n_stored);
    }
    std::iota(ordering.begin(), ordering.begin() + n_stored, 0);

    // Removing solo reads beforehand.
    size_t infront=0;
    for (size_t a=0; a<n_stored; ++a) {
        const auto& curlen=(remaining[a]=storage.get_len(a));
        if (curlen!=1) {
            continue;
        }
        
        if (static_cast<size_t>(*storage.get_start(a))!=a) {
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

        const auto curstart=storage.get_start(*maxIt);
        const auto curend=curstart + storage.get_len(*maxIt);
        std::swap(*maxIt, ordering[right-1]);
        --right;

        // Adding neighbors if they have not already been used somewhere else.
        size_t n_cluster=0;        
        for (auto curIt=curstart; curIt<curend; ++curIt) { 
            const auto& neighbor=*curIt;
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
