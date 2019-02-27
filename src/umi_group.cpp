#include "sarlacc.h"

#include "utils.h"
#include "DNA_input.h"
#include "sorted_trie.h"
#include "cluster_umis.h"

#include <vector>
#include <stdexcept>
#include <memory>
#include <deque>
#include <algorithm>

SEXP umi_group(SEXP umi1, SEXP thresh1, SEXP umi2, SEXP thresh2, SEXP pregroup) {
    BEGIN_RCPP

    // Checking the UMI values.
    auto seqs1=process_DNA_input(umi1);
    const size_t nseqs=seqs1->size();
    const int limit1=check_integer_scalar(thresh1, "threshold 1");

    const bool use_two=(umi2!=R_NilValue);
    std::unique_ptr<DNA_input> seqs2(nullptr);
    size_t nseq2=0;
    if (use_two) {
        seqs2=process_DNA_input(umi2);
        if (nseqs!=seqs2->size()) {
            throw std::runtime_error("'umi1' and 'umi2' should have the same length");
        }
    }
    const int limit2=check_integer_scalar(thresh2, "threshold 2");

    // Running through each group to define the current set. 
    Rcpp::List Pregroup(pregroup), output(Pregroup.size());
    for (size_t g=0; g<Pregroup.size(); ++g) {
        Rcpp::IntegerVector curgroup=Pregroup[g];
        const size_t curN=curgroup.size();

        if (curN==1) {
            output[g]=Rcpp::List::create(curgroup);
            continue;
        }

        std::vector<const char*> allseqs(curN);
        std::vector<int> alllens(curN), order(curN);
        value_store storage(curN);

        // Processing UMI1.
        seqs1->clear();
        for (size_t s=0; s<curN; ++s) {
            auto curseq=seqs1->get_persistent(curgroup[s] - 1);
            allseqs[s]=curseq.first;
            alllens[s]=curseq.second;
        }

        sorted_trie trie1(curN, allseqs.data(), alllens.data());
        sorted_trie::order(curN, allseqs.data(), alllens.data(), order.data());

        if (!use_two) {
            for (size_t s=0; s<curN; ++s) {
                auto o=order[s];
                storage[o]=trie1.find(allseqs[o], alllens[o], limit1);
            }

        } else { 
            value_store match_arrays(curN);

            // Saving matches for UMI1.
            for (size_t s=0; s<curN; ++s) {
                auto o=order[s];
                match_arrays[o]=trie1.find(allseqs[o], alllens[o], limit1);
                std::sort(match_arrays[o].begin(), match_arrays[o].end());
            }

            // Processing UMI2.
            seqs2->clear();
            for (size_t s=0; s<curN; ++s) {
                auto curseq=seqs2->get_persistent(curgroup[s] - 1);
                allseqs[s]=curseq.first;
                alllens[s]=curseq.second;
            }
            sorted_trie trie2(curN, allseqs.data(), alllens.data());
            sorted_trie::order(curN, allseqs.data(), alllens.data(), order.data());

            for (size_t s=0; s<curN; ++s) {
                auto o=order[s];
                const auto& matches=trie2.find(allseqs[o], alllens[o], limit2);
                std::deque<int> workspace;

                // Identifying the intersection for each read.
                const auto left=match_arrays[o].begin();
                const auto right=match_arrays[o].end();

                for (auto m2 : matches) {
                    auto it=std::lower_bound(left, right, m2);
                    if (it!=right && *it==m2) {
                        workspace.push_back(m2);
                    }
                }

                storage[o]=std::move(workspace);
            }
        }

        Rcpp::List curout=cluster_umis(storage);
        for (size_t i=0; i<curout.size(); ++i) {
            Rcpp::IntegerVector curvec=curout[i];
            for (auto& x : curvec) { x=curgroup[x]; }
            curout[i]=curvec;
        }
        output[g]=curout;
    }

    return output;
    END_RCPP
}


