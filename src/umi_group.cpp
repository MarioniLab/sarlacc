#include "sarlacc.h"
#include "utils.h"
#include "DNA_input.h"

struct trie_node {
    std::deque<int>* indices;
    std::vector<trie_node>* children;
    std::vector<int>* scores;
    size_t history;
    
    trie_node() : indices(NULL), children(NULL), scores(NULL), history(0) {}
    ~trie_node() {
        if (indices) { delete indices; }
        if (children) { delete children; }
        if (scores) { delete scores; }
    }

    bool dead_end () const {
        return (!indices && !children);
    }

    void _insert (const int index, const std::vector<char>& current, size_t position) {
        if (position==current.size()) {
            if (!indices) {            
                indices=new std::deque<int>();
            }
            indices->push_back(index);
            return;
        }

        if (!children) {
            children = new std::vector<trie_node>(NBASES);
        }

        switch (current[position]) {
            case 'A': 
                (*children)[0]._insert(index, current, position+1);
                break;
            case 'C': 
                (*children)[1]._insert(index, current, position+1);
                break;
            case 'G': 
                (*children)[2]._insert(index, current, position+1);
                break;
            case 'T':
                (*children)[3]._insert(index, current, position+1);
                break;
            case 'N':
                (*children)[4]._insert(index, current, position+1);
                break;
        }
        return;
    }

    // Top-level insert function. 
    void insert (int index, const std::vector<char>& current) {
        _insert(index, current, 0);
        return;
    }
    
    // Debugging function to explore tree shape.
    void _dump(int depth, char base) {
        for (int d=0; d<depth; ++d) {
            Rprintf("    ");
        }
        Rprintf("%c", base);
        if (indices) {
            Rprintf(": ");
            for (const auto& i : *indices) {
                Rprintf("%i ", i);
            }
        }
        Rprintf("\n");
        
        ++depth;
        if (children) {
            for (size_t i=0; i<NBASES; ++i) { 
                (*children)[i]._dump(depth, BASES[i]);
            }
        }
        return;
    }
    
    void dump() {
        _dump(0, 'x');        
    }

    // Computing the progressive levenshtein distance across the trie.
    void _find_within(std::deque<int>& collected, const std::vector<char>& current, size_t offset, 
                      char this_base, std::vector<int>& previous, int limit, size_t iter) {

        const size_t cur_len=current.size();
        if (!scores) {
            scores = new std::vector<int>(cur_len+1);
            *(scores->begin()) = *(previous.begin()) + 1;
            offset = 0;
        } else {
            if (scores->size() <= cur_len) {
                scores->resize(cur_len + 1);
            }
            if (history+1!=iter) {
                offset=0;
            }
        }
        history = iter;

        auto cur_space_it = scores->begin();
        auto last_space_it = previous.begin();
//      Rprintf("Last: ");
//        for (size_t i=0; i<=cur_len; ++i) {
//            Rprintf("%i, ", *(last_space_it+i));
//        }
//      Rprintf("\n");
//
//      Rprintf("This base is %c\n", this_base);

        // Filling out the latest row of the levenshtein distance matrix, picking up from the last memory (offset).
        // If the current base is 'N', this _always_ triggers a mismatch.
        if (this_base!='N') { 
            for (size_t i=offset+1; i<=cur_len; ++i) {
                *(cur_space_it+i) = std::min({ *(last_space_it + i) + 1, 
                                               *(cur_space_it + i - 1) + 1, 
                                               *(last_space_it + i - 1) + (current[i-1] == this_base ? 0 : 1) });
            }
        } else {
            for (size_t i=offset+1; i<=cur_len; ++i) {
                *(cur_space_it+i) = std::min({ *(last_space_it + i) + 1, 
                                               *(cur_space_it + i - 1) + 1, 
                                               *(last_space_it + i - 1) + 1 });
            }
        }

//        Rprintf("Current: ");
//        for (size_t i=0; i<=cur_len; ++i) {
//            Rprintf("%i, ", *(cur_space_it+i));
//        }
//        Rprintf("\n");
//        
//        for (size_t i=0; i<cur_len; ++i) {
//            Rprintf("%c", current[i]);
//        }
//        Rprintf("\n");

        // Adding indices to the result list, if they exist.
        const int cur_score = *(cur_space_it + cur_len); 
//        Rprintf("\tcurrent score is %i %i\n", cur_score, limit);        
        if (limit >= cur_score && indices) {  
            for (const auto& s : *indices) {
                collected.push_back(s);
            }
        }

        /* Checking if we should recurse through the children. 
         * This is not done if none of the scores across the final row of the
         * programming matrix are below or equal to limit, meaning that even 
         * a perfect match from here on in would always exceed limit. The 
         * current score doesn't count unless it's less than the limit,
         * as any extra extension would result in a +1 from that entry.
         */
        if (children) { 
            bool recurse_child = (limit > cur_score);
            if (!recurse_child && cur_len) { 
                for (auto rit = scores -> rbegin() + 1; rit!= (scores -> rend()); ++rit) {  
                    // Backwards, as scores at the end are probably lower and will trigger a break sooner.
                    if (*rit <= limit) {
                        recurse_child=true;
                        break;
                    }
                }
            }
            if (recurse_child) {
                for (size_t i=0; i<NBASES; ++i) {
                    auto& child=(*children)[i];
                    if (!child.dead_end()) { 
                        child._find_within(collected, current, offset, BASES[i], *scores, limit, iter);
                    }
                }
            }
        }
        return;
    }

    void find_within(std::deque<int>& collected, const std::vector<char>& current, int offset, int limit, size_t iter) {
        if (!scores) {
            scores = new std::vector<int>(current.size()+1);
        } else {
            scores->resize(current.size() + 1);
        }
        std::iota(scores->begin(), scores->end(), 0); // Initializing the first level of the levenshtein array here.

        if (children) {
            for (size_t i=0; i<NBASES; ++i) { 
                auto& child=(*children)[i];
                if (!child.dead_end()) { 
                    child._find_within(collected, current, offset, BASES[i], *scores, limit, iter);
                }
            }
        }

        return;
    }
};

struct sorted_trie {
    sorted_trie(SEXP s) : seqs(process_DNA_input(s)) {
        auto bcg_nspace=Rcpp::Environment::namespace_env("BiocGenerics");
        Rcpp::Function ordfun=bcg_nspace["order"];

        Rcpp::IntegerVector Order=ordfun(s);
        ordering.resize(Order.size());
        auto oIt=ordering.begin();
        for (auto o : Order) {
            (*oIt)=o-1; // zero-based indexing.
            ++oIt;
        }

        std::vector<char> current;
        current.reserve(50); // probably enough.
        
        for (auto o : ordering) {
            seqs->choose(o);
            const char * cur_str=seqs->cstring();
            const size_t cur_len=seqs->length();

            current.resize(cur_len);
            for (size_t i=0; i<cur_len; ++i) { 
                current[i]=seqs->decode(cur_str[i]);
            }

            toplevel.insert(o, current);
        }
        return;
    }

    Rcpp::List find_within(int limit) { 
        const size_t nseq=seqs->size();
        std::vector<char> current;
        current.reserve(50); // should probably be enough.

        Rcpp::List output(nseq);
        std::deque<int> collected;
        size_t counter=0;
        size_t last_o=0;

        for (auto o : ordering) {
            seqs->choose(o);
            const char * cur_str=seqs->cstring();
            const size_t cur_len=seqs->length();

            // Figuring out the common prefix.
            size_t idx=0;
            while (idx < current.size() && idx < cur_len && seqs->decode(cur_str[idx])==current[idx]) {
                ++idx;
            }
            const size_t common=idx;

            // If it's fully common, this implies that this string is the same as the previous string, 
            // so we use the previous results and skip to avoid redundant work.
            if (common==cur_len && cur_len==current.size() && counter) { 
                output[o]=output[last_o];
                continue;
            }

            // Replacing the remaining elements in 'current'.
            current.resize(cur_len);
            while (idx < cur_len) {
                current[idx]=seqs->decode(cur_str[idx]);
                ++idx;                
            }

//            for (size_t i=0; i<cur_len; ++i) {
//                Rprintf("%c", current[i]);
//            }
//            Rprintf("\n");
//            Rprintf("Common is %i\n", common);

            // Searching through the trie.
            toplevel.find_within(collected, current, common, limit, counter);

            // Storing the results (getting back to 1-based indexing).
            for (auto& x : collected) { ++x; }
            output[o]=Rcpp::IntegerVector(collected.begin(), collected.end());
            collected.clear();

            ++counter;
            last_o=o;
        }
        
        return output;
    }

    void dump() {
        toplevel.dump();
    }

private:
    trie_node toplevel;
    std::vector<size_t> ordering;
    std::unique_ptr<DNA_input> seqs;
};

/***************************************
 * R-level functions start here.
 ***************************************/

SEXP umi_group (SEXP umi, SEXP threshold) {
    BEGIN_RCPP
    int limit=check_integer_scalar(threshold, "distance threshold");
    sorted_trie dic(umi);
    return dic.find_within(limit);
    END_RCPP
}
