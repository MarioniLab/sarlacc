#include "sarlacc.h"
#include "utils.h"
#include "DNA_input.h"
#include "sorted_trie.h"
#include "umi_clusterer.h"
#include "value_store.h"

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
    const int limit2=check_integer_scalar(thresh1, "threshold 2");

    // Defining persistent and reusable arrays for intermediate variable-length constructs.
    std::vector<const char*> allseqs;
    std::vector<int> alllens, order;
    umi_clusterer clusterer;
    
    value_store<int, std::deque<int> > match_arrays;
    std::vector<int> workspace;

    // Running through each group to define the current set. 
    Rcpp::List Pregroup(pregroup), output(Pregroup.size());
    for (size_t g=0; g<Pregroup.size(); ++g) {
        Rcpp::IntegerVector curgroup=Pregroup[g];
        const size_t curN=curgroup.size();

        if (curN==1) {
            output[g]=curgroup;
            continue;
        }

        if (curN > allseqs.size()) {
            allseqs.resize(curN);
            alllens.resize(curN);
            order.resize(curN);
        }

        // Processing UMI1.
        for (size_t s=0; s<curN; ++s) {
            auto curseq=seqs1->get(curgroup[s]);
            allseqs[s]=curseq.first;
            alllens[s]=curseq.second;
        }
        sorted_trie trie1(curN, allseqs.data(), alllens.data());
        sorted_trie::order(curN, allseqs.data(), alllens.data(), order.data());

        if (!use_two) {
            for (auto o : order) {
                auto matches=trie1.find(allseqs[o], alllens[o], limit1);
                clusterer.storage.add(matches.begin(), matches.end(), o);
            }

        } else { 
            // Saving matches for UMI1.
            size_t used=0;
            for (auto o : order) {
                auto matches=trie1.find(allseqs[o], alllens[o], limit1);
                match_arrays.add(matches.begin(), matches.end(), o);
                auto startIt=match_arrays.get_start_unsafe(o);
                std::sort(startIt, startIt + match_arrays.get_len(o));
            }

            // Processing UMI2.
            for (size_t s=0; s<curN; ++s) {
                auto curseq=seqs2->get(curgroup[s]);
                allseqs[s]=curseq.first;
                alllens[s]=curseq.second;
            }
            sorted_trie trie2(curN, allseqs.data(), alllens.data());
            sorted_trie::order(curN, allseqs.data(), alllens.data(), order.data());

            for (auto o : order) {
                auto matches=trie2.find(allseqs[o], alllens[o], limit2);
                const size_t n_match_one=match_arrays.get_len(o);

                const size_t max_size=std::min(n_match_one, matches.size());
                if (max_size > workspace.size()) {
                    workspace.resize(max_size);
                }

                // Identifying the intersection for each read.
                const auto left=match_arrays.get_start(o);
                const auto right=left + n_match_one;
                int counter=0;

                for (auto m2 : matches) {
                    auto it=std::lower_bound(left, right, m2);
                    if (it!=right && *it==m2) {
                        workspace[counter]=m2;
                        ++counter;
                    }
                }

                clusterer.storage.add(workspace.begin(), workspace.begin() + counter);
            }

            match_arrays.clear();
        }

        output[g]=clusterer.cluster();
        clusterer.storage.clear();
    }

    return output;
    END_RCPP
}


