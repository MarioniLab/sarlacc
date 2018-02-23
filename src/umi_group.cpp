#include "sarlacc.h"

const char START='x';

struct trie_node {
    std::deque<int>* indices;
    std::vector<trie_node>* children;
    
    trie_node() : indices(NULL), children(NULL) {}
    ~trie_node() {
        if (indices) { delete indices; }
        if (children) { delete children; }
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
            children = new std::vector<trie_node>(4);
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
            (*children)[0]._dump(depth, 'A');
            (*children)[1]._dump(depth, 'C');
            (*children)[2]._dump(depth, 'G');
            (*children)[3]._dump(depth, 'T');
        }
        return;
    }
    
    void dump() {
        _dump(0, START);        
    }

    // Computing the progressive levenshtein distance across the trie.
    void _find_within(std::deque<int>& collected, const std::vector<char>& current, char this_base, 
                      std::deque<std::vector<int> >& running, int limit) const {

        if (!indices && !children) {
            // Don't bother computing if I have no indices or children.
            return; 
        }
        
//      Rprintf("Searching for %c\n", this_base);        
        const size_t cur_len=current.size();
        auto last_space_it = running.back().begin();
        running.push_back(std::vector<int>(cur_len+1));
        auto cur_space_it=running.back().begin();

//      Rprintf("Last: ");
//        for (size_t i=0; i<=cur_len; ++i) {
//            Rprintf("%i, ", *(last_space_it+i));
//        }
//      Rprintf("\n");

//      Rprintf("This base is %c\n", this_base);
        // Picking up from the last memory.
        *(cur_space_it)=*(last_space_it)+1;
        for (size_t i=1; i<=cur_len; ++i) {
            *(cur_space_it+i) = std::min({ *(last_space_it + i) + 1, 
                                           *(cur_space_it + i - 1) + 1, 
                                           *(last_space_it + i - 1) + (current[i-1] == this_base ? 0 : 1) });
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
            if (!recurse_child) { 
                for (size_t i=0; i<cur_len; ++i) {
                    if (*(cur_space_it + i) <= limit) {
                        recurse_child=true;
                        break;
                    }
                }
            }
            if (recurse_child) {
                (*children)[0]._find_within(collected, current, 'A', running, limit);
                (*children)[1]._find_within(collected, current, 'C', running, limit);
                (*children)[2]._find_within(collected, current, 'G', running, limit);
                (*children)[3]._find_within(collected, current, 'T', running, limit);
            }
        }

        // Popping off the vector we added before returning to the calling function.
        running.pop_back();
        return;
    }

    void find_within(std::deque<int>& collected, const std::vector<char>& current, int limit) const {
        std::deque<std::vector<int> > running(1, std::vector<int>(current.size()+1));
        std::iota(running.back().begin(), running.back().end(), 0); // Initializing the first level of the levenshtein array here.

        if (children) {
            (*children)[0]._find_within(collected, current, 'A', running, limit);
            (*children)[1]._find_within(collected, current, 'C', running, limit);
            (*children)[2]._find_within(collected, current, 'G', running, limit);
            (*children)[3]._find_within(collected, current, 'T', running, limit);
        }

        return;
    }
};

struct sorted_trie {
    sorted_trie(const XStringSet_holder * p) : ptr(p) {
        const size_t nseq=get_length_from_XStringSet_holder(ptr);

        for (size_t s=0; s<nseq; ++s) {
            auto cur_data=get_elt_from_XStringSet_holder(ptr, s);
            const char * cur_str=cur_data.ptr;
            const size_t cur_len=cur_data.length;

            std::vector<char> current(cur_str, cur_str + cur_len);
            for (auto& x : current) { x=DNAdecode(x); }

            toplevel.insert(s, current);
        }
        return;
    }

    void find_within(std::deque<int>& results, int index, int limit) {
        auto cur_data=get_elt_from_XStringSet_holder(ptr, index);
        const char * cur_str=cur_data.ptr;
        const size_t cur_len=cur_data.length;
        
        std::vector<char> current(cur_str, cur_str + cur_len);
        for (auto& x : current) { x=DNAdecode(x); }

        toplevel.find_within(results, current, limit);
        return;
    }

    void dump() {
        toplevel.dump();
    }
    
    trie_node toplevel;
    const XStringSet_holder * ptr;
};

SEXP umi_group_x (SEXP umi, SEXP threshold) {
    BEGIN_RCPP
        
    // Checking inputs.       
    auto seq=hold_XStringSet(umi); 
    const size_t nseq=get_length_from_XStringSet_holder(&seq);
    int limit=check_integer_scalar(threshold, "distance threshold");

    sorted_trie dic(&seq);
    std::deque<int> collected;
    Rcpp::List output(nseq);

    for (size_t i=0; i<nseq; ++i) {
        dic.find_within(collected, i, limit);
        for (auto& x : collected) { ++x; }
        output[i]=Rcpp::IntegerVector(collected.begin(), collected.end());
        collected.clear();
    }
   
    return output;
    END_RCPP
}
