#include "sarlacc.h"
#include "DNA_input.h"
#include "utils.h"
#include "masker.h"

#include <seqan/align.h>
#include <seqan/graph_msa.h>

SEXP quick_msa(SEXP sequences, SEXP groupings, SEXP threshold, SEXP qualities, SEXP submat, SEXP open, SEXP ext) {
    BEGIN_RCPP
    
    // Checking inputs.
    auto all_seq=process_DNA_input(sequences);
    const size_t nseq=all_seq->size();

    Rcpp::List Groups(groupings);
    const size_t ngroups=Groups.size();

    // Setting up the scoring scheme.
    seqan::Score<double, seqan::ScoreMatrix<seqan::Dna5, seqan::Default> > score_scheme(
        check_numeric_scalar(ext, "gap extension penalty"),
        check_numeric_scalar(open, "gap opening penalty")
    );

    Rcpp::NumericMatrix Substitution(submat);
    if (Substitution.nrow() != Substitution.ncol()) { 
        throw std::runtime_error("substitution matrix should be square");
    }

    const size_t nvalues=seqan::ValueSize<seqan::Dna5>::VALUE;
    for (size_t i = 0; i < nvalues; ++i) {
        for (size_t j = 0; j < nvalues; ++j) {
            setScore(score_scheme, seqan::Dna5(i), seqan::Dna5(j), Substitution(i,j));
        }
    }

    // Checking whether we want to mask by qualities.
    const double max_error=check_numeric_scalar(threshold, "masking threshold");
    bool use_quals=(qualities!=R_NilValue && !ISNA(max_error));
    Rcpp::List Quals;
    if (use_quals) {
        Quals=Rcpp::List(qualities);
        if (Quals.size()!=nseq){
            throw std::runtime_error("length of quality list is not the same as number of sequences");
        }
    }

    const size_t maxseqwidth=get_max_width(all_seq.get());
    masker seqM(max_error, use_quals ? maxseqwidth : 0);
    unmasker seqU(use_quals ? maxseqwidth : 0);

    // Iterating through all groups of reads to perform the MSA.
    Rcpp::List output(ngroups);
    for (size_t g=0; g<ngroups; ++g) {
        Rcpp::IntegerVector curgroup(Groups[g]);
        const size_t cursize=curgroup.size();
       
        // Skipping or breaking if there aren't enough reads in the group. 
        if (cursize==0L) {
            throw std::runtime_error("empty group of reads for multiple sequence alignment");
        } else if (cursize==1L) {
            auto curpair=all_seq->get(curgroup[0]-1);
            output[g]=Rcpp::StringVector::create(Rcpp::String(curpair.first));
            continue;
        }

        // Filling up the alignment object, possibly with masked sequences.
        seqan::Align<seqan::String<seqan::Dna5> > align;
        seqan::resize(seqan::rows(align), cursize);
            
        all_seq->clear();
        if (use_quals) {
            for (size_t i=0; i<cursize; ++i) {
                auto curpair=all_seq->get(curgroup[i]-1);
                Rcpp::NumericVector curqual(Quals[curgroup[i]-1]);
                Rcpp::String out=seqM.mask(curpair.first, curpair.second, curqual.begin());
                seqan::assignSource(seqan::row(align, i), out.get_cstring());
            }
        } else {
            for (size_t i=0; i<cursize; ++i) {
                auto curpair=all_seq->get(curgroup[i]-1);
                seqan::assignSource(seqan::row(align, i), curpair.first);
            }
        }

        // Running the MSA (T-coffee).
        seqan::globalMsaAlignment(align, score_scheme);

        // Storing the results (possibly with unmasking).
        Rcpp::StringVector curout(curgroup.size());
        for (size_t i=0; i<cursize; ++i) {
            std::stringstream out;
            out << seqan::row(align, i);

            if (use_quals) {
                auto tmp=out.str();
                auto original=all_seq->get(curgroup[i]-1);
                curout[i]=seqU.unmask(tmp.c_str(), tmp.size(), original.first, original.second);
            } else {
                curout[i]=out.str();
            }
        }
        output[g]=curout;
    }

    return output;
    END_RCPP
}
