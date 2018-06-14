#include "sarlacc.h"
#include "DNA_input.h"
#include "utils.h"

// Constants.
const double max_error=0.99999999, min_error=0.00000001;

Rcpp::String errorsToString (size_t len, const std::vector<double>& errorprobs, std::vector<char>& working) {
    if (working.size() <= len) {
        working.resize(len+1);
    }
    for (size_t i=0; i<len; ++i) {
        const int to_ascii=std::min(
            std::round(-10 * errorprobs[i] / std::log(10)),
            93.0 /* corresponding to a quality of '~' */ );

        working[i]=char(to_ascii + 33);
    }

    working[len]='\0';
    return Rcpp::String(working.data());
}

struct data_holder {
    std::vector<char> buffer;
    std::vector<double> scores;
    std::vector<int> incidences;
    std::vector<double> errorprobs;

    void expand(size_t alignwidth) {
        if (incidences.size() <= alignwidth) {
            incidences.resize(alignwidth);
            scores.resize(alignwidth * NBASES);
            buffer.resize(alignwidth + 1);
            errorprobs.resize(alignwidth);
        }

        std::fill(incidences.begin(), incidences.begin() + alignwidth, 0);
        std::fill(scores.begin(), scores.begin() + alignwidth*NBASES, 0);
        std::fill(buffer.begin(), buffer.begin() + alignwidth + 1, '\0');
        std::fill(errorprobs.begin(), errorprobs.begin() + alignwidth, 0);
        return;
    }
};

/* Function to compute consensus without any quality scores.
 * I've split it into an internal function, which performs the calculation per alignment;
 * and an external function, which accepts a list of alignments for efficient looping.
 */

size_t internal_create_consensus_basic(SEXP alignments, const double mincov, const double pseudo_denom, data_holder& storage) {
    auto all_aln=process_DNA_input(alignments);
    const size_t naligns=all_aln->size();
    const double pseudo_num=pseudo_denom/NBASES;

    const size_t alignwidth=check_alignment_width(all_aln.get());
    storage.expand(alignwidth);    

    // Counting the number of occurrences of each base at each position.
    for (size_t a=0; a<naligns; ++a) {
        all_aln->choose(a);
        const char* aln_str=all_aln->cstring();
        auto sIt=storage.scores.begin();

        for (size_t i=0; i<alignwidth; ++i) {
            const char curbase=all_aln->decode(aln_str[i]);
            switch (curbase) { 
                case 'A': case 'a':
                    ++(*sIt);
                    break;
                case 'C': case 'c':
                    ++(*(sIt+1));
                    break;
                case 'G': case 'g':
                    ++(*(sIt+2));
                    break;
                case 'T': case 't':
                    ++(*(sIt+3));
                    break;
            }
            
            // We use a separate vector for holding incidences, just in case
            // there are "N"'s (such that the sum of ACTG counts != incidences).
            if (curbase!='-') {
                ++(storage.incidences[i]);
            }
            sIt+=NBASES;
        }
    }

    // Constructing the consensus sequence.
    auto sIt=storage.scores.begin();
    auto cIt=storage.buffer.begin();
    auto eIt=storage.errorprobs.begin();

    for (size_t i=0; i<alignwidth; ++i, sIt+=NBASES) {
        if (storage.incidences[i] < double(naligns)*mincov) {
             continue;
        }

        // Picking the most abundant base.
        auto maxed=std::max_element(sIt, sIt+NBASES);
        (*cIt)=BASES[maxed - sIt];
        ++cIt;

        // Denominator computed from the non-N bases only; N's neither help nor hinder, as the data is missing.
        const double total=std::accumulate(sIt, sIt+NBASES, 0);
        const double correct_prob=(*maxed + pseudo_num)/(total + pseudo_denom);
        *eIt=std::log1p(-correct_prob); // remember, _error_ probabilities.
        ++eIt;
    }

    return cIt - storage.buffer.begin();
}

SEXP create_consensus_basic(SEXP alignments, SEXP min_cov, SEXP pseudo_count) { 
    BEGIN_RCPP
    const double mincov=check_numeric_scalar(min_cov, "minimum coverage");
    const double pseudo_denom=check_numeric_scalar(pseudo_count, "pseudo count");

    data_holder storage;
    size_t conlen=internal_create_consensus_basic(alignments, mincov, pseudo_denom, storage);

    auto eIt=storage.errorprobs.begin();
    return Rcpp::List::create(Rcpp::String(storage.buffer.data()), Rcpp::NumericVector(eIt, eIt+conlen));
    END_RCPP
}

SEXP create_consensus_basic_loop(SEXP alignments, SEXP min_cov, SEXP pseudo_count) { // loop over a list of alignments.
    BEGIN_RCPP
    Rcpp::List aln_list(alignments);
    const double mincov=check_numeric_scalar(min_cov, "minimum coverage");
    const double pseudo_denom=check_numeric_scalar(pseudo_count, "pseudo count");

    const size_t nalign=aln_list.size();
    data_holder storage;
    Rcpp::StringVector all_cons(nalign);
    Rcpp::StringVector all_qual(nalign);
    std::vector<char> qual_buffer;
    
    for (size_t i=0; i<nalign; ++i) {
        size_t conlen=internal_create_consensus_basic(aln_list[i], mincov, pseudo_denom, storage);
        all_cons[i]=Rcpp::String(storage.buffer.data());
        all_qual[i]=errorsToString(conlen, storage.errorprobs, qual_buffer);
    }

    return Rcpp::List::create(all_cons, all_qual);
    END_RCPP
}

/* This creates a consensus where the choice of base is aware of the different qualities across different reads. 
 * Similarly, the output quality score calculation will take that into account.
 * Again, I've split it into an internal function, which performs the calculation per alignment;
 * and an external function, which accepts a list of alignments for efficient looping.
 */

size_t internal_create_consensus_quality(SEXP alignments, const double mincov, SEXP qualities, data_holder& storage) {
    auto all_aln=process_DNA_input(alignments);
    const size_t naligns=all_aln->size();
    const size_t alignwidth=check_alignment_width(all_aln.get());
    storage.expand(alignwidth);

    Rcpp::List qual(qualities); // need the numeric qualities as they are decoded differently depending on whether they are Solexa or Phred.
    const size_t nquals=qual.size();
    if (nquals!=naligns) {
        throw std::runtime_error("alignments and qualities have different numbers of entries");
    }

    // Running through each entry.
    for (size_t a=0; a<naligns; ++a) {
        all_aln->choose(a);
        const char* astr=all_aln->cstring();

        Rcpp::NumericVector curqual(qual[a]);
        auto sIt=storage.scores.begin();
        int position=0;

        for (size_t i=0; i<alignwidth; ++i, sIt+=NBASES) { // leave sIt here to ensure it runs even when 'continue's.
            const char curbase=all_aln->decode(astr[i]);
            if (curbase=='-') {
                continue;
            }
            ++(storage.incidences[i]);

            if (position >= curqual.size()) {
                throw std::runtime_error("quality vector is shorter than the alignment sequence");
            }
            double newqual=curqual[position];
            if (newqual>max_error) {
                newqual=max_error;
            } else if (newqual<min_error) {
                newqual=min_error;
            }

            const double right=std::log1p(-newqual); // using base e to make life easier.
            const double wrong=std::log(newqual/(NBASES-1));
            ++position;

            // Adding the log-probabilities for each base. This is safe with 'N's, as all bases get +=wrong.
            for (size_t b=0; b<NBASES; ++b) {
                *(sIt+b) += (curbase==BASES[b] ? right : wrong);
            }
        }

        if (position!=curqual.size()) { 
            throw std::runtime_error("quality vector is longer than the alignment sequence");
        }
    }

    // Constructing the consensus sequence.
    auto sIt=storage.scores.begin();
    auto cIt=storage.buffer.begin();
    auto eIt=storage.errorprobs.begin();

    for (size_t i=0; i<alignwidth; ++i, sIt+=NBASES) { // leave sIt here to ensure it runs even when 'continue's are triggered.
        if (storage.incidences[i] < double(naligns)*mincov) {
             continue;
        }

        // Choosing the base with the highest probability.
        auto maxed=std::max_element(sIt, sIt+NBASES);
        (*cIt)=BASES[maxed - sIt];
        ++cIt;
        
        // Summing probabilities to get the denominator of the probabilities.
        std::sort(sIt, sIt+NBASES);
        double denom=*sIt, error;
        for (size_t b=1; b<NBASES; ++b) {
            const double leftover=*(sIt+b) - denom;
            denom+=R::log1pexp(leftover);

            if (b==NBASES-2) { // i.e., excluding the maximum log-probability at position 'upto-1'.
                error=denom;
            }
        }

        (*eIt)=error - denom;
        ++eIt;
    }

    return cIt - storage.buffer.begin();
}

SEXP create_consensus_quality(SEXP alignments, SEXP min_cov, SEXP qualities) { // loop over a list of alignments.
    BEGIN_RCPP
    const double mincov=check_numeric_scalar(min_cov, "minimum coverage");

    data_holder storage;
    size_t conlen=internal_create_consensus_quality(alignments, mincov, qualities, storage);

    auto eIt=storage.errorprobs.begin();
    return Rcpp::List::create(Rcpp::String(storage.buffer.data()), Rcpp::NumericVector(eIt, eIt+conlen));
    END_RCPP
}

SEXP create_consensus_quality_loop(SEXP alignments, SEXP min_cov, SEXP qualities) { // loop over a list of alignments.
    BEGIN_RCPP
    Rcpp::List aln_list(alignments);
    Rcpp::List qual_list(qualities);
    const double mincov=check_numeric_scalar(min_cov, "minimum coverage");

    const size_t nalign=aln_list.size();
    data_holder storage;
    Rcpp::StringVector all_cons(nalign);
    Rcpp::StringVector all_qual(nalign);
    std::vector<char> qual_buffer;
    
    for (size_t i=0; i<nalign; ++i) {
        size_t conlen=internal_create_consensus_quality(aln_list[i], mincov, qual_list[i], storage);
        all_cons[i]=Rcpp::String(storage.buffer.data());
        all_qual[i]=errorsToString(conlen, storage.errorprobs, qual_buffer);
    }
    
    return Rcpp::List::create(all_cons, all_qual);
    END_RCPP
}
