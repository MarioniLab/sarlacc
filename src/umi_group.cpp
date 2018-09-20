#include "sarlacc.h"
#include "utils.h"
#include "DNA_input.h"
#include "sorted_trie.h"
#include "umi_clusterer.h"

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
    std::vector<int> match_from_one(1000), workspace;
    std::vector<size_t> match_starts, match_num; 
    std::vector<const char*> allseqs;
    std::vector<int> alllens, order;
    umi_clusterer clusterer;

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
                clusterer.add(matches.begin(), matches.end());
            }

        } else { 
            // Saving matches for UMI1.
            if (curN > match_starts.size()) {
                match_starts.resize(curN);
                match_num.resize(curN);
            }

            size_t used=0;
            for (auto o : order) {
                auto matches=trie1.find(allseqs[o], alllens[o], limit1);
                match_starts[o]=used;
                match_num[o]=matches.size();

                const size_t required=used+matches.size();
                if (required > match_from_one.size()) {
                    match_from_one.resize(required);
                }
                std::copy(matches.begin(), matches.end(), match_from_one.begin() + used);
                std::sort(match_from_one.begin() + used, match_from_one.begin() + required);
                used=required;
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
                const auto left=match_from_one.begin() + match_starts[o];
                const auto right=left + match_num[o];

                const size_t max_size=std::min(match_num[o], matches.size());
                if (max_size > workspace.size()) {
                    workspace.resize(max_size);
                }

                // Identifying the intersection for each read.
                int counter=0;
                for (auto i : matches) {
                    auto it=std::lower_bound(left, right, i);
                    if (it!=right && *it==i) {
                        workspace[counter]=i;
                        ++counter;
                    }
                }

                clusterer.add(workspace.begin(), workspace.begin() + counter);
            }            
        }

        Rcpp::List clusters=clusterer.cluster();
        for (size_t i=0; i<clusters.size(); ++i) {
            Rcpp::IntegerVector current=clusters[i];
            for (auto& val : current) {
                val=order[val];
            }
            clusters[i]=current;
        }
        
        output[g]=clusters;
        clusterer.clear();
    }

    return output;
    END_RCPP
}


