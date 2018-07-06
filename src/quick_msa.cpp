#include "sarlacc.h"
#include "DNA_input.h"
#include <seqan/align.h>
#include <seqan/graph_msa.h>

SEXP quick_msa(SEXP sequences, SEXP groupings) {
    BEGIN_RCPP
    
    // Checking inputs.
    auto all_seq=process_DNA_input(sequences);
    const size_t nseq=all_seq->size();

    Rcpp::List Groups(groupings);
    const size_t ngroups=Groups.size();

    Rcpp::List output(ngroups);
    for (size_t g=0; g<ngroups; ++g) {
        Rcpp::IntegerVector curgroup(Groups[g]);
        const size_t cursize=curgroup.size();
        all_seq->clear();

        seqan::Align<seqan::String<seqan::Dna> > align;
        seqan::resize(seqan::rows(align), cursize);
        for (size_t i=0; i<cursize; ++i) {
            auto curpair=all_seq->get(curgroup[i]-1);
            seqan::assignSource(seqan::row(align, i), curpair.first);
        }
        
        seqan::globalMsaAlignment(align, seqan::SimpleScore());

        Rcpp::StringVector curout(curgroup.size());
        for (size_t i=0; i<cursize; ++i) {
            std::stringstream out;
            out << seqan::row(align, i);
            curout[i]=out.str();
        }
        output[g]=curout;
    }

    return output;
    END_RCPP
}
