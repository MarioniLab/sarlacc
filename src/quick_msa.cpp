#include "sarlacc.h"

#include "DNA_input.h"
#include "utils.h"

#include <seqan/align.h>
#include <seqan/graph_msa.h>

#include <sstream>

typedef seqan::Dna5String TSource;
typedef seqan::Score<int, seqan::Simple> TScore;
typedef seqan::StringSet<TSource, seqan::Dependent<> > TStringSet;

SEXP quick_msa(SEXP groupings, SEXP sequences, SEXP match, SEXP mismatch, SEXP gapExtension, SEXP gapOpening, SEXP bandwidth) {
    BEGIN_RCPP
    
    // Checking inputs.
    auto all_seq=process_DNA_input(sequences);
    const size_t nseq=all_seq->size();
    Rcpp::List Groups(groupings);
    const size_t ngroups=Groups.size();

    // Setting up SeqAn options.
    seqan::MsaOptions<TSource, TScore> msaOpt;
    msaOpt.sc = TScore(
        check_numeric_scalar(match, "match score"), 
        check_numeric_scalar(mismatch, "mismatch score"), 
        check_numeric_scalar(gapExtension, "gap extension score"), 
        check_numeric_scalar(gapOpening, "gap opening score")
    );
    msaOpt.isDefaultPairwiseAlignment = false;
    msaOpt.method = 0; // global alignment
    msaOpt.pairwiseAlignmentMethod = 2; // banded
    msaOpt.bandWidth = check_integer_scalar(bandwidth, "bandwidth");

    Rcpp::List output(ngroups);

    for (size_t g=0; g<ngroups; ++g) {
        Rcpp::IntegerVector curgroup(Groups[g]);
        const size_t cursize=curgroup.size();

        // Skipping groups with only one read.
        if (cursize==0) {
            output[g]=Rcpp::StringVector(0);
        } else if (cursize==1) {
            const char* seqstr=all_seq->get(curgroup[0]-1).first;
            output[g]=Rcpp::StringVector::create(Rcpp::String(seqstr));
            continue;
        }

        // Proceeding to T-coffee.
        seqan::Align<TSource> align;
        seqan::resize(seqan::rows(align), cursize);
        for (size_t i=0; i<cursize; ++i) {
            auto curpair=all_seq->get(curgroup[i]-1);
            seqan::assignSource(seqan::row(align, i), curpair.first);
        }

        // Copied out of "graph_align_tcoffee_msa.h" to enable passing of alignment options.
        TStringSet sequenceSet = seqan::stringSet(align);
        seqan::Graph<seqan::Alignment<TStringSet, void, seqan::WithoutEdgeId> > gAlign(sequenceSet);

        seqan::String<seqan::String<char> > sequenceNames;
        seqan::resize(sequenceNames, seqan::length(sequenceSet), seqan::String<char>("tmpName"));
        seqan::globalMsaAlignment(gAlign, sequenceSet, sequenceNames, msaOpt);
        seqan::convertAlignment(gAlign, align);

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
