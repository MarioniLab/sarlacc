#include "cluster_umis.h"

#include <vector>
#include <algorithm>
#include <stdexcept>

Rcpp::List cluster_umis(const value_store& storage) {
    const size_t n_stored=storage.size();

	std::vector<size_t> ordering(n_stored), remaining(n_stored);
    std::deque<int> per_cluster_values;
    std::deque<Rcpp::IntegerVector> output;
    
    if (n_stored > ordering.size()) { 
        ordering.resize(n_stored);
        remaining.resize(n_stored);
    }
    std::iota(ordering.begin(), ordering.begin() + n_stored, 0);

    // Removing solo reads beforehand.
    size_t infront=0, n_out=0;
    for (size_t a=0; a<n_stored; ++a) {
        const auto& curlen=(remaining[a]=storage[a].size());
        if (curlen > 1) {
            continue;
        } else if (curlen==1) {
            if (static_cast<size_t>(storage[a].front())!=a) {
                throw std::runtime_error("single-read groups should contain only the read itself");
            }

            Rcpp::IntegerVector tosave(1, a);
            if (n_out < output.size()) {
                output[n_out]=tosave;
            } else {
                output.push_back(tosave);
            }
            ++n_out;
        } else if (curlen==0) {
            throw std::runtime_error("zero length read group");
        }

        std::swap(ordering[a], ordering[infront]);
        ++infront;
    }

    // Keeping the top nodes, accounting for updated counts after processing prior nodes.
    auto left=ordering.begin() + infront, right=ordering.begin() + n_stored;
    while (left < right) {

        // Wiping out empty nodes with an effective pop_front().
        for (auto it=left; it!=right; ++it) {
            if (remaining[*left]==0) {
                std::swap(*left, *it);
                ++left;
            }
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

        auto curstart=storage[*maxIt].begin();
        auto curend=storage[*maxIt].end();
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
            for (auto& next : storage[neighbor]) {
                auto& nextremain=remaining[next];
                if (nextremain > 0) { 
                    --nextremain;
                }
            }
        }

        Rcpp::IntegerVector tosave(per_cluster_values.begin(), per_cluster_values.begin() + n_cluster);
        if (n_out < output.size()) {
            output[n_out]=tosave;
        } else {
            output.push_back(tosave);
        }
        ++n_out;
    }

    return Rcpp::List(output.begin(), output.begin() + n_out);
}


