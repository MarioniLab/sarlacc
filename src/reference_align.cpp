#include "reference_align.h"

reference_align::reference_align (size_t reflen, const char * refseq, Rcpp::NumericVector qualities, double go, double ge) : 
        rlen(reflen), rseq(refseq), gap_open(go+ge), gap_ext(ge), 
        dpmatrix(nrows*(rlen+1)), affine_left(nrows), 
        backtrack_start(rlen+1), backtrack_end(rlen+1) 
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
    size_t newsize=nrows*(rlen+1);
    if (newsize > dpmatrix.size()) { 
        dpmatrix.resize(newsize);
        affine_left.resize(nrows);
    }

    // Filling the first column of the DP matrix.
    auto location=dpmatrix.begin();
    if (local) {
        // If local, filling with zeroes (no penalty for vertical gaps at start).
        std::fill(location, location+nrows, std::make_pair(up, 0));
    } else {
        // If global, filling more conventionally with vertical gap penalties.
        location->first=up;
        location->second=0;
        auto loccopy=location+1;
        for (size_t i=1; i<nrows; ++i, ++loccopy) {
            loccopy->first=up;
            loccopy->second=- gap_open - gap_ext * (i-1);
        }
    }
    
    // Filling the initial horizontal affine penalties.
    std::fill(affine_left.begin(), affine_left.begin() + nrows, R_NegInf);

    for (size_t col=1; col<rlen; ++col) {
        location+=nrows;
        align_column(location, affine_left.begin(), rseq[col-1], len, seq, qual, false);
    }
        
    if (rlen) {
        location+=nrows;
        // If 'local', then last=true to avoid vertical gap penalties.
        align_column(location, affine_left.begin(), rseq[rlen-1], len, seq, qual, local);
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
    backtracked=false;
    return (location+len)->second;
}

void reference_align::align_column(std::deque<dpentry>::iterator storage, std::deque<double>::iterator last_horizontal, char reference, 
        size_t len, const char* seq, const char* qual, bool last) 
{
    auto lastcol=storage - len - 1; // Assume that we're past the first column of the DP matrix.
    auto lagging=storage;
    auto lagging_past=lastcol;

    // Compute stats for the first row of the DP matrix separately.
    storage->first=left;
    storage->second = lastcol->second - (lastcol->first == left ? gap_ext : gap_open);
    ++storage;
    ++lastcol;
    ++last_horizontal;

    double vert_gap_open=(last ? 0 : gap_open);
    double vert_gap_ext=(last ? 0 : gap_ext);
    double last_vertical=R_NegInf;

    for (size_t i=0; i<len; ++i, ++storage, ++lastcol, ++lagging, ++lagging_past, ++last_horizontal) {
        // Horizontal gap opening cost. We need to compute the cost
        // of extending an already opened gap that was suboptimal  
        // in the last column.
        double horiz_gap = lastcol->second - (lastcol->first==left ? gap_ext : gap_open); 
        double previous_horiz_gap = *last_horizontal - gap_ext;
        bool reopen_horiz = previous_horiz_gap > horiz_gap;
        if (reopen_horiz) {
            horiz_gap = previous_horiz_gap;
        }
        *last_horizontal = horiz_gap;

        // Vertical gap opening cost. Again, we need to compute the
        // cost of extending an already opened gap that was suboptimal
        // in the last row.
        double vert_gap = lagging->second - (lagging->first==up ? vert_gap_ext : vert_gap_open);
        double previous_vert_gap = last_vertical - vert_gap_ext;
        bool reopen_vert = previous_vert_gap > vert_gap;
        if (reopen_vert) {
            vert_gap = previous_vert_gap;
        }
        last_vertical = vert_gap;

        // (Mis)match cost.
        double match = lagging_past->second + compute_cost(reference, seq[i], qual[i]);
        
        if (match > horiz_gap && match > vert_gap) {
            storage->first=diag;
            storage->second=match;
        } else if (horiz_gap > vert_gap) {
            storage->first=left;
            storage->second=horiz_gap;
            if (reopen_horiz) { // Alter history so that horizontal gap is already opened.
                lastcol->first=left;
            }
        } else {
            storage->first=up;
            storage->second=vert_gap;
            if (reopen_vert) { // Again, alter history so that vertical gap is already opened.
                lagging->first=up;
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

void reference_align::backtrack(bool include_gaps) {
    if (!aligned) {
        throw std::runtime_error("cannot backtrack without alignment");
    }
    backtracked=true;

    auto location=dpmatrix.begin() + nrows * rlen; // start at the start of the last column.
    if (include_gaps) {
        backtrack_end[rlen]=nrows;
    } else {
        // Deciding whether or not to include end gaps.
        size_t currow=nrows - 1;
        auto loccopy=location + currow;
        while (currow > 0 && loccopy->first==up) {
            --loccopy;
            --currow;
        }
        backtrack_end[rlen]=currow+1;
    }

    for (size_t i=rlen; i>0; --i) { // starting from the last column and working forwards. 
        size_t currow=backtrack_end[i]-1;
        auto loccopy=location + currow;
        while (currow > 0 && loccopy->first==up) {
            --loccopy;
            --currow;
        }

        location-=nrows;
        backtrack_start[i]=currow;
        if (loccopy->first==left) {
            backtrack_end[i-1]=currow+1;
        } else {
            backtrack_end[i-1]=currow;
        }
    }
   
    if (include_gaps) {
        backtrack_start[0]=0;
    } else {
        backtrack_start[0]=backtrack_end[0]-1;
    }

#ifdef DEBUG
    for (size_t i=0; i<= rlen; ++i) {
        Rprintf("%i %i %i\n", i, backtrack_start[i], backtrack_end[i]);
    }
#endif    
    return;
}

std::pair<size_t, size_t> reference_align::map (size_t ref_start, size_t ref_end) const {
    if (!backtracked) {
        throw std::runtime_error("cannot extract regions without backtracking");
    }

    // Protect against empty references.
    if (!rlen) {
        return std::make_pair(0, 0);
    }

    // +1 to convert to DP matrix column index; -1 to get back to sequence index from DP row index.
    return std::make_pair(
        std::max(size_t(1), backtrack_start[ref_start+1]) - 1,
        std::max(size_t(1), backtrack_end[ref_end]) - 1 // no +1, as we want the column corresponding to the before-the-end base.
    );
}
