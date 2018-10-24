#include "reference_align.h"

reference_align::reference_align (size_t reflen, const char * refseq, Rcpp::NumericVector qualities, double go, double ge) : 
        rlen(reflen), rseq(refseq), gap_open(go+ge), gap_ext(ge), 
        dpmatrix(nrows*(rlen+1)), 
        backtrack_start(rlen+1), backtrack_end(rlen+1) 
{
    create_qualities(qualities);
    return;
};

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

double reference_align::align(size_t len, const char* seq, const char* qual) {
    nrows=len+1;
    size_t newsize=nrows*(rlen+1);
    if (newsize > dpmatrix.size()) { 
        dpmatrix.resize(newsize);
    }

    // Filling the first column of the DP matrix with zeroes (no penalty for vertical gap openings at start).
    auto location=dpmatrix.begin();
    std::fill(location, location+nrows, std::make_pair(up, 0));

    for (size_t col=1; col<rlen; ++col) {
        location+=nrows;
        align_column(location, rseq[col-1], len, seq, qual, false);
    }
        
    location+=nrows;
    if (rlen) {
        align_column(location, rseq[rlen-1], len, seq, qual, true);
    }

    aligned=true;
    backtracked=false;
    return (location+len)->second;
}

void reference_align::align_column(std::deque<dpentry>::iterator storage, char reference, size_t len, const char* seq, const char* qual, bool last) {
    auto lastcol=storage - len - 1; // Assume that we're past the first column of the DP matrix.
    auto lagging=storage;
    auto lagging_past=lastcol;

    // Compute stats for the first row of the DP matrix separately.
    storage->first=left;
    storage->second = lastcol->second - (lastcol->first == left ? gap_ext : gap_open);
    ++storage;
    ++lastcol;

    double vert_gap_open=(last ? 0 : gap_open);
    double vert_gap_ext=(last ? 0 : gap_ext);

    for (size_t i=0; i<len; ++i, ++storage, ++lastcol, ++lagging, ++lagging_past) {
        // Horizontal gap opening cost.
        double horiz_gap = lastcol->second - (lastcol->first == left ? gap_ext : gap_open);

        // Vertical gap opening cost.
        double vert_gap = lagging->second - (lagging->first == up ? vert_gap_ext : vert_gap_open);

        // (Mis)match cost.
        double match = lagging_past->second + compute_cost(reference, seq[i], qual[i]);

        if (match > horiz_gap && match > vert_gap) {
            storage->first=diag;
            storage->second=match;
        } else if (horiz_gap > vert_gap) {
            storage->first=left;
            storage->second=horiz_gap;
        } else {
            storage->first=up;
            storage->second=vert_gap;
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
