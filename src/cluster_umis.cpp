#include "sarlacc.h"

SEXP cluster_umis(SEXP groupings) {
    BEGIN_RCPP

    Rcpp::List Groupings(groupings);
    const size_t nreads=Groupings.size();

    std::deque<size_t> ordering(nreads);
    std::iota(ordering.begin(), ordering.end(), 0);
    
    std::vector<size_t> remaining(nreads);
    std::vector<Rcpp::IntegerVector> storage(nreads);
    for (size_t i=0; i<nreads; ++i) {
        storage[i]=Groupings[i];
        remaining[i]=storage[i].size();
    }

    // Removing solo reads beforehand.
    std::deque<Rcpp::IntegerVector> output;
    {
        size_t infront=0;
        for (size_t o=0; o<ordering.size(); ++o) {
            size_t& curdex=ordering[o];
            const size_t& current=remaining[curdex];
            if (current!=1) {
                continue;
            }

            if (storage[curdex][0]-1!=curdex) {
                throw std::runtime_error("single-read groups should contain only the read itself");
            }
            output.push_back(storage[curdex]);
            std::swap(curdex, ordering[infront]);
            ++infront;
        }

        for (size_t i=0; i<infront; ++i) {
            ordering.pop_front();
        }
    }
                
    // Keeping the top nodes, accounting for updated counts after processing prior nodes.
    std::deque<int> collected;
    while (ordering.size()) {
       
        // Wiping out empty nodes.
        size_t infront=0;
        for (auto& o : ordering) { 
            const size_t& current=remaining[o];
            if (current==0) {
                std::swap(o, ordering[infront]);
                ++infront;
            } 
        }
        for (size_t i=0; i<infront; ++i) {
            ordering.pop_front();
        }
        if (ordering.empty()) {
            break;
        }

        // Finding the node with the maximum size.
        auto& maxed=*std::max_element(ordering.begin(), ordering.end(), 
            [&] (const size_t& left, const size_t& right) -> bool {
                if (remaining[left]==remaining[right]) {
                    return left < right;
                }
                return remaining[left] < remaining[right];
            }
        );
        std::swap(maxed, ordering.back());
        const auto& curgroup=storage[ordering.back()];
        ordering.pop_back();
        collected.clear();
       
        // Adding neighbors if they have not already been used somewhere else.
        for (auto neighbor : curgroup) {
            --neighbor; // zero indexing.
            auto& remain=remaining[neighbor];
            if (remain==0) {
                continue;
            }
            collected.push_back(neighbor+1);
            remain=0;

            // Also decrementing the counts of all neighbors of neighbors.
            const auto& neighgroup=storage[neighbor];
            for (auto nextneighbor : neighgroup) {
                auto& nextremain=remaining[nextneighbor-1];
                if (nextremain > 0) { 
                    --nextremain;
                }
            }
        }

        output.push_back(Rcpp::IntegerVector(collected.begin(), collected.end()));
    }
    
    return Rcpp::List(output.begin(), output.end());
    END_RCPP
}
