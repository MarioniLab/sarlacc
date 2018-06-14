#include "sarlacc.h"
#include "DNA_input.h"

// Constants.
const std::vector<char> BASES={'A', 'C', 'G', 'T'};
const int NBASES=BASES.size();
const double max_error=0.99999999, min_error=0.00000001;

/* Function to compute consensus without any quality scores.
 * I've split it into an internal function, which performs the calculation per alignment;
 * and an external function, which accepts a list of alignments for efficient looping.
 */

std::pair<Rcpp::String, Rcpp::String> internal_create_consensus_basic(SEXP alignments, const double mincov, const double pseudo_denom) {
    auto all_aln=process_DNA_input(alignments);
    const size_t naligns=all_aln->size();
    const size_t alignwidth=check_alignment_width(all_aln.get());
    const double pseudo_num=pseudo_denom/NBASES;
    
    std::vector<int> scores(alignwidth*NBASES), incidences(alignwidth);

    // Counting the number of occurrences of each base at each position.
    for (size_t a=0; a<naligns; ++a) {
        all_aln->choose(a);
        const char* aln_str=all_aln->cstring();
        auto sIt=scores.begin();

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
                ++incidences[i];
            }
            sIt+=NBASES;
        }
    }

    // Constructing the consensus sequence.
    std::vector<char> consensus(alignwidth+1, '\0');
    std::vector<char> qual_str(alignwidth+1, '\0');
    auto sIt=scores.begin();
    auto cIt=consensus.begin();
    auto qIt=qual_str.begin();

    for (size_t i=0; i<alignwidth; ++i, sIt+=NBASES) {
        if (incidences[i] < double(naligns)*mincov) {
             continue;
        }

        // Picking the most abundant base.
        auto maxed=std::max_element(sIt, sIt+NBASES);
        (*cIt)=BASES[maxed - sIt];
        ++cIt;

        // Denominator computed from the non-N bases only; N's neither help nor hinder, as the data is missing.
        const double total=std::accumulate(sIt, sIt+NBASES, 0);
        const double correct_prob=(*maxed + pseudo_num)/(total + pseudo_denom);
        const double logerr=std::log1p(-correct_prob);
        (*qIt) = char(-10 * logerr/std::log(10) + 33);
        ++qIt;
    }

    return std::make_pair(Rcpp::String(consensus.data()), Rcpp::String(qual_str.data()));
}

SEXP create_consensus_basic(SEXP alignments, SEXP min_cov, SEXP pseudo_count) { // loop over a list of alignments.
    BEGIN_RCPP
    Rcpp::List aln_list(alignments);
    const double mincov=check_numeric_scalar(min_cov, "minimum coverage");
    const double pseudo_denom=check_numeric_scalar(pseudo_count, "pseudo count");

    const size_t nalign=aln_list.size();
    Rcpp::StringVector all_cons(nalign);
    Rcpp::StringVector all_seqs(nalign);
    
    for (size_t i=0; i<nalign; ++i) {
        auto out=internal_create_consensus_basic(aln_list[i], mincov, pseudo_denom);
        all_cons[i]=out.first;
        all_seqs[i]=out.second;
    }

    return Rcpp::List::create(all_cons, all_seqs);
    END_RCPP
}

/* This creates a consensus where the choice of base is aware of the different qualities across different reads. 
 * Similarly, the output quality score calculation will take that into account.
 * Again, I've split it into an internal function, which performs the calculation per alignment;
 * and an external function, which accepts a list of alignments for efficient looping.
 */

std::pair<Rcpp::String, Rcpp::String> internal_create_consensus_quality(SEXP alignments, const double mincov, SEXP qualities) {
    auto all_aln=process_DNA_input(alignments);
    const size_t naligns=all_aln->size();
    const size_t alignwidth=check_alignment_width(all_aln.get());

    Rcpp::List qual(qualities); // need the numeric qualities as they are decoded differently depending on whether they are Solexa or Phred.
    const size_t nquals=qual.size();
    if (nquals!=naligns) {
        throw std::runtime_error("alignments and qualities have different numbers of entries");
    }

    std::vector<int> incidences(alignwidth);
    std::vector<double> scores(alignwidth*NBASES);

    // Running through each entry.
    for (size_t a=0; a<naligns; ++a) {
        all_aln->choose(a);
        const char* astr=all_aln->cstring();

        Rcpp::NumericVector curqual(qual[a]);
        auto sIt=scores.begin();
        int position=0;

        for (size_t i=0; i<alignwidth; ++i, sIt+=NBASES) { // leave sIt here to ensure it runs even when 'continue's.
            const char curbase=all_aln->decode(astr[i]);
            if (curbase=='-') {
                continue;
            }
            ++incidences[i];

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
    std::vector<char> consensus(alignwidth+1, '\0');
    std::vector<char> qual_str(alignwidth+1, '\0');
    auto sIt=scores.begin();
    auto cIt=consensus.begin();
    auto qIt=qual_str.begin();

    for (size_t i=0; i<alignwidth; ++i, sIt+=NBASES) { // leave sIt here to ensure it runs even when 'continue's are triggered.
        if (incidences[i] < double(naligns)*mincov) {
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

        (*qIt)=char((error - denom) / std::log(10) * -10 + 33);
        ++qIt;
    }

    return std::make_pair(Rcpp::String(consensus.data()), Rcpp::String(qual_str.data()));
}

SEXP create_consensus_quality(SEXP alignments, SEXP min_cov, SEXP qualities) { // loop over a list of alignments.
    BEGIN_RCPP
    Rcpp::List aln_list(alignments);
    Rcpp::List qual_list(qualities);
    const double mincov=check_numeric_scalar(min_cov, "minimum coverage");

    const size_t nalign=aln_list.size();
    Rcpp::StringVector all_cons(nalign);
    Rcpp::StringVector all_seqs(nalign);
    
    for (size_t i=0; i<nalign; ++i) {
        auto out=internal_create_consensus_quality(aln_list[i], mincov, qual_list[i]);
        all_cons[i]=out.first;
        all_seqs[i]=out.second;
    }

    return Rcpp::List::create(all_cons, all_seqs);
    END_RCPP
}
