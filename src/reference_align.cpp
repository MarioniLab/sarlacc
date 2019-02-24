#include "reference_align.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

reference_align::reference_align (size_t reflen, const char * refseq, Rcpp::NumericVector qualities, double go, double ge) : 
        rlen(reflen), rseq(refseq), gap_open(go+ge), gap_ext(ge), 
        scores(nrows), left_jump_scores(nrows), left_jump_points(nrows), directions(nrows*(rlen+1))
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
        left_jump_scores.resize(nrows);
        left_jump_points.resize(nrows);
        directions.resize(nrows * (rlen + 1));
    }

    // Filling the first column of the DP matrix (negative values mean 'up' direction).
    std::fill(directions.begin(), directions.begin()+nrows, -1);
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
    std::fill(left_jump_scores.begin(), left_jump_scores.begin() + nrows, R_NegInf);
    std::fill(left_jump_points.begin(), left_jump_points.begin() + nrows, 0);

    // Processing all columns corresponding to reference positions.
    auto last_dir=directions.begin(); 
    for (size_t col=1; col<rlen; ++col) {
        align_column(last_dir, col-1, len, seq, qual, false);
        last_dir+=nrows;
    }
        
    if (rlen) {
        // If 'local', then last=true to avoid vertical gap penalties.
        align_column(last_dir, rlen-1, len, seq, qual, local);
    }

#ifdef DEBUG
    for (size_t i=0;i<nrows; ++i) {
        auto loc=directions.begin() + i;
        for(size_t j=0;j<=rlen; ++j) {
            Rprintf("%i\t", *loc);
            loc+=nrows;
        }
        Rprintf("\n");
    }
#endif

    aligned=true;
    return scores[len];
}

void reference_align::align_column(std::deque<int>::iterator last_direction,
        size_t pos, size_t len, const char* seq, const char* qual, bool last) 
{
    const char reference=rseq[pos];
    auto current_direction=last_direction+nrows;

    // Compute stats for the first row of the DP matrix separately.
    // Note that 'scores', upon input, holds the scores of the *last* column;
    // but upon output, it holds the scores of the *current* column.
    double lagging_last = scores[0];
    scores[0] -= (*last_direction > 0 ? gap_ext : gap_open);
    *current_direction=1;

    double vert_gap_open=(last ? 0 : gap_open);
    double vert_gap_ext=(last ? 0 : gap_ext);
    double up_jump_score=R_NegInf;
    size_t up_jump_point=0;

    for (size_t i=1; i<=len; ++i) {
        // Horizontal gap opening cost. We need to compute the cost
        // of extending an already opened gap that was suboptimal  
        // in the previous columns (i.e., along the current row).
        double horiz_gap = scores[i] - (*(last_direction+i) > 0 ? gap_ext : gap_open);
        
        double& previous_horiz_gap = left_jump_scores[i];
        previous_horiz_gap -= gap_ext;
        size_t left_step=1;
        if (previous_horiz_gap > horiz_gap) {
            left_step=1 + pos - left_jump_points[i];
            horiz_gap=previous_horiz_gap;
        } else {
            previous_horiz_gap=horiz_gap;
            left_jump_points[i]=pos;
        }

        // Vertical gap opening cost. Again, we need to compute the
        // cost of extending an already opened gap that was suboptimal
        // in the previous rows (i.e., along the current column).
        double vert_gap = scores[i-1] - (*(current_direction+i-1) < 0 ? vert_gap_ext : vert_gap_open);

        up_jump_score -= vert_gap_ext;
        size_t up_step=1;
        if (up_jump_score > vert_gap) {
            up_step=1 + i - up_jump_point;
            vert_gap=up_jump_score;
        } else {
            up_jump_score=vert_gap;
            up_jump_point=i;
        }

        // (Mis)match cost. Note 'i' is +1 from the base position of the new sequence, hence the use of 'i-1'.
        // We also update the lagging score from the last column.
        double match = lagging_last + compute_cost(reference, seq[i-1], qual[i-1]);
        lagging_last = scores[i];

        auto& curdir=*(current_direction+i);
        auto& curscore=scores[i];
        if (match > horiz_gap && match > vert_gap) {
            curdir=0;
            curscore=match;
        } else if (horiz_gap > vert_gap) {
            curscore=horiz_gap;
            curdir=left_step;
        } else {
            curscore=vert_gap;
            curdir=up_step;
            curdir*=-1;
        }

//        Rprintf("\tchoices: %.5f %.5f %.5f\n", match, horiz_gap, vert_gap);
//        Rprintf("\t\t%i %.5f %.5f %i %.5f %i\n", curdir, curscore, left_jump_scores[i], left_jump_points[i], up_jump_score, up_jump_point);
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

/* Backtracking template function, to avoid having to rewrite
 * the same code in similar contexts. 
 */

template<class U>
void reference_align::backtrack(U& store) const {
    if (!aligned) {
        throw std::runtime_error("cannot backtrack without alignment");
    }

    auto location=directions.begin() + nrows * rlen; // start at the start of the last column.
    size_t currow=nrows-1;

    // Starting from the last column and working forwards.
    // Note that we don't fill i=0, but this doesn't correspond to a real base anyway.
    for (size_t i=rlen; i>0; --i) { 
        while (currow > 0) {
            auto curdir=*(location + currow);
            if (curdir >= 0) {
                break;
            }

            while (curdir < 0) {
                store.move_up(i, currow);
                --currow;
                ++curdir;
            }
        }

        auto curdir=*(location+currow);
        if (curdir==0) {
            store.move_diag(i, currow);
            --currow;
        } else {
            while (curdir > 0) {
                store.move_left(i, currow);
                --curdir;
            }
        }

        location-=nrows;
    }

    while (currow > 0) {
        store.move_up(0, currow);
        --currow;
    }
    return;
}

void reference_align::fill_map(querymap& path_store) const {
    struct mapper {
        mapper(size_t rlen, std::deque<size_t>& s, std::deque<size_t>& e) : backtrack_start(s), backtrack_end(e) {
            backtrack_start.resize(rlen+1);
            backtrack_end.resize(rlen+1);
        }
        void move_up(size_t i, size_t currow) {}
        void move_diag(size_t i, size_t currow) {
            // Only include a query position if it's an actual match.
            backtrack_start[i]=currow;
            backtrack_end[i]=currow+1;
            return;
        }
        void move_left(size_t i, size_t currow) {
            // A horizontal gap is not 'aligned' to the current reference base.
            backtrack_start[i]=currow+1;
            backtrack_end[i]=currow+1;
            return;            
        }
        std::deque<size_t>& backtrack_start;
        std::deque<size_t>& backtrack_end;
    };

    mapper store(rlen, path_store.starts, path_store.ends);
    backtrack(store);
   
#ifdef DEBUG
    for (size_t i=0; i<= rlen; ++i) {
        Rprintf("%i %i %i\n", i, backtrack_start[i], backtrack_end[i]);
    }
#endif    
    return;
}

std::pair<size_t, size_t> reference_align::querymap::operator()(size_t ref_start, size_t ref_end) {
    if (starts.size()!=ends.size()) {
        throw std::runtime_error("invalid input");
    }
    if (starts.size()<=1) {
        return std::make_pair(0, 0);
    }

    return std::make_pair(
        // +1 to convert to DP matrix column index; -1 to get back to sequence index from DP row index.
        starts[ref_start+1] - 1,
        // no +1, as we want the column corresponding to the before-the-end base.
        ends[ref_end] - 1 
    );
}

void reference_align::fill_strings(std::vector<char>& rwork, std::vector<char>& qwork, const char* qseq) const {
    struct stringer {
        stringer(std::vector<char>& r, std::vector<char>& q, const char* rs, const char* qs) : rwork(r), qwork(q), rseq(rs), qseq(qs) {
            // Vector resize() will not reduce capacity, avoid reallocs.
            rwork.resize(1, '\0');
            qwork.resize(1, '\0');
            return;
        }
        void move_up(size_t i, size_t currow) {
            rwork.push_back('-');
            qwork.push_back(qseq[currow-1]);
            return;
        }
        void move_left(size_t i, size_t currow) {
            rwork.push_back(rseq[i-1]);
            qwork.push_back('-');
            return;
        }
        void move_diag(size_t i, size_t currow) {
            rwork.push_back(rseq[i-1]);
            qwork.push_back(qseq[currow-1]);
            return;
        }
        std::vector<char>& rwork;
        std::vector<char>& qwork;
        const char* rseq, * qseq;
    };

    stringer store(rwork, qwork, rseq, qseq);
    backtrack(store);

    std::reverse(rwork.begin(), rwork.end());
    std::reverse(qwork.begin(), qwork.end());
    return;
}


