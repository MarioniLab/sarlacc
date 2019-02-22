#include "reference_align.h"

reference_align::reference_align (size_t reflen, const char * refseq, Rcpp::NumericVector qualities, double go, double ge) : 
        rlen(reflen), rseq(refseq), gap_open(go+ge), gap_ext(ge), 
        scores(nrows), affine_left(nrows), directions(nrows*(rlen+1))
{
    create_qualities(qualities);
    return;
}

/* Given the error probability 'epsilon' for the observed base, we have the score as:
 *
 * \log_2(gamma_{x,y} * (1 - epsilon) * n + (1 - gamma_{x,y}) * epsilon * (n/(n-1)))
 *
 * where n = 4. We don't have to worry about the error probability for the known reference.
 */
void reference_align::create_qualities (Rcpp::NumericVector qualities) {
    quality_encoding Q(qualities);
    offset=Q.lowest();
    auto errors=Q.get_errors();
    available=errors.size();

    precomputed_match.resize(4);
    precomputed_mismatch.resize(4);
    for (size_t i=0; i<4; ++i) {
        auto& mscores=precomputed_match[i];
        auto& mmscores=precomputed_mismatch[i];
        mscores.resize(available);
        mmscores.resize(available);

        double gamma_xy=1.0/(i+1.0);
        double gamma_xy_1m = 1 - gamma_xy;
        constexpr double n = 4;

        for (size_t j=0; j<available; ++j) {
            const double epsilon=errors[j];
            mscores[j]=std::log(gamma_xy * (1 - epsilon) * n + gamma_xy_1m * epsilon * (n/(n-1)))/M_LN2;
        }

        // gamma_xy => 1 - gamma_xy for a mismatch.
        for (size_t j=0; j<available; ++j) {
            const double epsilon=errors[j];
            mmscores[j]=std::log(gamma_xy_1m * (1 - epsilon) * n + gamma_xy * epsilon * (n/(n-1)))/M_LN2;
        }
    }
    
    return;
}

double reference_align::align(size_t len, const char* seq, const char* qual, bool local) {
    nrows=len+1;
    if (nrows > scores.size()) {
        scores.resize(nrows);
        affine_left.resize(nrows);
        directions.resize(nrows * (rlen + 1));
    }

    // Filling the first column of the DP matrix.
    std::fill(directions.begin(), directions.begin()+nrows, up);
    if (local) {
        // If local, filling with zeroes (no penalty for vertical gaps at start).
        std::fill(scores.begin(), scores.begin()+nrows, 0);
    } else {
        // If global, filling more conventionally with vertical gap penalties.
        scores[0]=0;
        for (size_t i=1; i<nrows; ++i) {
            scores[i]=- gap_open - gap_ext * (i-1);
        }
    }
    
    // Filling the initial horizontal affine penalties.
    std::fill(affine_left.begin(), affine_left.begin() + nrows, R_NegInf);

    // Processing all columns corresponding to reference positions.
    auto last_dir=directions.begin(); 
    for (size_t col=1; col<rlen; ++col) {
        align_column(last_dir, rseq[col-1], len, seq, qual, false);
        last_dir+=nrows;
    }
        
    if (rlen) {
        // If 'local', then last=true to avoid vertical gap penalties.
        align_column(last_dir, rseq[rlen-1], len, seq, qual, local);
    }

#ifdef DEBUG
    for (size_t i=0;i<nrows; ++i) {
        auto loc=dpmatrix.begin() + i;
        for(size_t j=0;j<=rlen; ++j) {
            Rprintf("%.3f\t", loc->second);
            loc+=nrows;
        }
        Rprintf("\n");
    }
#endif

    aligned=true;
    return scores[len];
}

void reference_align::align_column(std::deque<DIRECTION>::iterator last_direction,
        char reference, size_t len, const char* seq, const char* qual, bool last) 
{
    auto current_direction=last_direction+nrows;

    // Compute stats for the first row of the DP matrix separately.
    // Note that 'scores', upon input, holds the scores of the *last* column;
    // but upon output, it holds the scores of the *current* column.
    double lagging_last = scores[0];
    scores[0] -= (*last_direction == left ? gap_ext : gap_open);
    *current_direction=left;

    double vert_gap_open=(last ? 0 : gap_open);
    double vert_gap_ext=(last ? 0 : gap_ext);
    double affine_up=R_NegInf;

    for (size_t i=1; i<=len; ++i) {
        // Horizontal gap opening cost. We need to compute the cost
        // of extending an already opened gap that was suboptimal  
        // in the last column.
        double horiz_gap = scores[i] - (*(last_direction+i)==left ? gap_ext : gap_open); 
        double previous_horiz_gap = affine_left[i] - gap_ext;
        bool reopen_horiz = previous_horiz_gap > horiz_gap;
        if (reopen_horiz) {
            horiz_gap = previous_horiz_gap;
        }
        affine_left[i] = horiz_gap;

        // Vertical gap opening cost. Again, we need to compute the
        // cost of extending an already opened gap that was suboptimal
        // in the last row.
        auto prev=i-1;
        double vert_gap = scores[prev] - (*(current_direction+prev)==up ? vert_gap_ext : vert_gap_open);
        double previous_vert_gap = affine_up - vert_gap_ext;
        bool reopen_vert = previous_vert_gap > vert_gap;
        if (reopen_vert) {
            vert_gap = previous_vert_gap;
        }
        affine_up = vert_gap;

        // (Mis)match cost. Note 'i' is +1 from the base position of the new sequence, hence the use of 'prev'.
        // We also update the lagging score from the last column.
        double match = lagging_last + compute_cost(reference, seq[prev], qual[prev]);
        lagging_last = scores[i];

        auto& curdir=*(current_direction+i);
        auto& curscore=scores[i];
        if (match > horiz_gap && match > vert_gap) {
            curdir=diag;
            curscore=match;
        } else if (horiz_gap > vert_gap) {
            curdir=left;
            curscore=horiz_gap;
            if (reopen_horiz) { // Alter history so that horizontal gap is already opened.
                *(last_direction+i)=left;
            }
        } else {
            curdir=up;
            curscore=vert_gap;
            if (reopen_vert) { // Again, alter history so that vertical gap is already opened.
                *(current_direction+prev)=up;
            }
        }
    }

    return;
}
   
// We assume that the observed sequence can only be 'A', 'C', 'G' or 'T'.
double reference_align::compute_cost (char ref, char obs, char qual) const {
    switch (ref) {
        case 'A': case 'C': case 'G': case 'T':
            return precomputed_cost(1, (ref==obs), qual);
        case 'M':
            return precomputed_cost(2, (ref=='A' || ref=='C'), qual);
        case 'R': 
            return precomputed_cost(2, (ref=='A' || ref=='G'), qual);
        case 'W':
            return precomputed_cost(2, (ref=='A' || ref=='T'), qual);
        case 'S':
            return precomputed_cost(2, (ref=='C' || ref=='G'), qual);
        case 'Y':
            return precomputed_cost(2, (ref=='C' || ref=='T'), qual);
        case 'K':
            return precomputed_cost(2, (ref=='G' || ref=='T'), qual);
        case 'V':
            return precomputed_cost(3, (ref!='T'), qual);
        case 'H':
            return precomputed_cost(3, (ref!='G'), qual);
        case 'D':
            return precomputed_cost(3, (ref!='C'), qual);
        case 'B':
            return precomputed_cost(3, (ref!='A'), qual);
        case 'N':
            return precomputed_cost(4, true, qual);
    }
    throw std::runtime_error("unrecognized base in reference sequence");
}

double reference_align::precomputed_cost (int mode, bool matched, char qual) const {
    if (qual < offset) {
        throw std::runtime_error("quality cannot be lower than smallest encoded value");
    }
    size_t location=qual - offset;
    if (location >= available) {
        location=available - 1;
    }
    
    auto& origin=(matched ? precomputed_match : precomputed_mismatch);
    return origin[mode - 1][location];
}

void reference_align::backtrack(std::deque<size_t>& backtrack_start, std::deque<size_t>& backtrack_end) const {
    if (!aligned) {
        throw std::runtime_error("cannot backtrack without alignment");
    }
    backtrack_start.resize(rlen+1);
    backtrack_end.resize(rlen+1);

    auto location=directions.begin() + nrows * rlen; // start at the start of the last column.
    size_t currow=nrows-1;

    // Starting from the last column and working forwards.
    // Note that we don't fill i=0, but this doesn't correspond to a real base anyway.
    for (size_t i=rlen; i>0; --i) { 
        while (currow > 0 && *(location+currow)==up) {
            --currow;
        }

        // Only include a query position if it's an actual match;
        // a horizontal gap is not 'aligned' to the current reference base.
        if (*(location+currow)==left) { 
            backtrack_start[i]=currow+1;
            backtrack_end[i]=currow+1;
        } else {
            backtrack_start[i]=currow;
            backtrack_end[i]=currow+1;
            --currow;
        }
        
        location-=nrows;
    }
   
#ifdef DEBUG
    for (size_t i=0; i<= rlen; ++i) {
        Rprintf("%i %i %i\n", i, backtrack_start[i], backtrack_end[i]);
    }
#endif    
    return;
}

void reference_align::backtrack(std::deque<char>& ref_str, std::deque<char>& query_str, const char* qseq) const {
    if (!aligned) {
        throw std::runtime_error("cannot backtrack without alignment");
    }
    ref_str.clear();
    query_str.clear();

    auto location=directions.begin() + nrows * rlen; // start at the start of the last column.
    size_t currow=nrows-1;

    // Starting from the last column and working forwards.
    // Note that we don't fill i=0, but this doesn't correspond to a real base anyway.
    for (size_t i=rlen; i>0; --i) { 
        while (currow > 0 && *(location+currow)==up) {
            ref_str.push_front('-');
            query_str.push_front(qseq[currow-1]);
            --currow;
        }

        // Only include a query position if it's an actual match;
        // a horizontal gap is not 'aligned' to the current reference base.
        if (*(location+currow)==left) {
            ref_str.push_front(rseq[i-1]);
            query_str.push_front('-');
        } else {
            ref_str.push_front(rseq[i-1]);
            query_str.push_front(qseq[currow-1]);
            --currow;
        }
        
        location-=nrows;
    }
   
    return;
}

std::pair<size_t, size_t> reference_align::map (const std::deque<size_t>& backtrack_start, const std::deque<size_t>& backtrack_end, size_t ref_start, size_t ref_end) {
    if (backtrack_start.size()!=backtrack_end.size()) {
        throw std::runtime_error("invalid input");
    }
    if (backtrack_start.size()<=1) {
        return std::make_pair(0, 0);
    }

    return std::make_pair(
        // +1 to convert to DP matrix column index; -1 to get back to sequence index from DP row index.
        backtrack_start[ref_start+1] - 1,
        // no +1, as we want the column corresponding to the before-the-end base.
        backtrack_end[ref_end] - 1 
    );
}
