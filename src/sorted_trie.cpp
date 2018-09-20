#include "sorted_trie.h"

const std::vector<char> BASES={'A', 'C', 'G', 'T', 'N'};
const int NBASES=BASES.size();

const int MULT=2;
const int EDIT=MULT;
int get_edit_score (char left, char right) {
    // If the current base is 'N', this _always_ adds half a mismatch.
    // This reduces the penalty of a full mismatch without awarding ambiguity.
    if (left=='N' || right=='N') { return 1; }
    if (left==right) { return 0; }
    return EDIT;
}

/*****************************************
 * Internal class for nodes of the trie. *
 *****************************************/

sorted_trie::trie_node::trie_node() : indices(NULL), children(NULL), scores(NULL), history(0) {}

sorted_trie::trie_node::~trie_node() {
	if (indices) { delete indices; }
	if (children) { delete children; }
	if (scores) { delete scores; }
}

bool sorted_trie::trie_node::dead_end () const {
	return (!indices && !children);
}

void sorted_trie::trie_node::insert (const int index, const std::vector<char>& current, size_t position) {
	if (position==current.size()) {
		if (!indices) {
			indices=new std::deque<int>(1, index);
		} else {
    		indices->push_back(index);
        }
		return;
	}

	if (!children) {
		children = new std::vector<trie_node>(NBASES);
	}

	switch (current[position]) {
		case 'A':
			(*children)[0].insert(index, current, position+1);
			break;
		case 'C':
			(*children)[1].insert(index, current, position+1);
			break;
		case 'G':
			(*children)[2].insert(index, current, position+1);
			break;
		case 'T':
			(*children)[3].insert(index, current, position+1);
			break;
		case 'N':
			(*children)[4].insert(index, current, position+1);
			break;
	}
	return;
}

// Top-level insert function.
void sorted_trie::trie_node::insert (int index, const std::vector<char>& current) {
	insert(index, current, 0);
	return;
}

// Debugging function to explore tree shape.
void sorted_trie::trie_node::dump(int depth, char base) {
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
			(*children)[i].dump(depth, BASES[i]);
		}
	}
	return;
}

void sorted_trie::trie_node::dump() {
	dump(0, 'x');
}

// Computing the progressive levenshtein distance across the trie.
void sorted_trie::trie_node::find_within(std::deque<int>& collected, const std::vector<char>& current, size_t offset, char this_base, std::vector<int>& previous, int limit, size_t iter) {
	const size_t cur_len=current.size();
	if (!scores) {
		scores = new std::vector<int>(cur_len+1);
		*(scores->begin()) = *(previous.begin()) + EDIT;
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
//    Rprintf("Last: ");
//      for (size_t i=0; i<=cur_len; ++i) {
//          Rprintf("%i, ", *(last_space_it+i));
//      }
//    Rprintf("\n");
//
//    Rprintf("This base is %c\n", this_base);

    // Filling out the latest row of the levenshtein distance matrix, picking up from the last memory (offset).
    for (size_t i=offset+1; i<=cur_len; ++i) {
        *(cur_space_it+i) = std::min({ *(last_space_it + i) + EDIT,
                                       *(cur_space_it + i - 1) + EDIT,
                                       *(last_space_it + i - 1) + get_edit_score(current[i-1], this_base) });
    }

//    Rprintf("Current: ");
//    for (size_t i=0; i<=cur_len; ++i) {
//        Rprintf("%i, ", *(cur_space_it+i));
//    }
//    Rprintf("\n");
//
//    for (size_t i=0; i<cur_len; ++i) {
//        Rprintf("%c", current[i]);
//    }
//    Rprintf("\n");

    // Adding indices to the result list, if they exist.
    const int cur_score = *(cur_space_it + cur_len);
//    Rprintf("\tcurrent score is %i %i\n", cur_score, limit);
    if (limit >= cur_score && indices) {
        for (const auto& s : *indices) {
            collected.push_back(s);
        }
    }

    /* Checking if we should recurse through the children.
     * This is not done if all of the scores across the final row of the programming matrix are above 'limit',
     * meaning that even a perfect match from here on in would always exceed 'limit'.
     * The current position (bottom-right on the array) only counts if it's lless than 'limit - EDIT',
     * as any extra extension would result in a '+EDIT' to the score from that entry.
     */
    if (children) {
        bool recurse_child = (limit >= cur_score + EDIT);
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
                    child.find_within(collected, current, offset, BASES[i], *scores, limit, iter);
                }
            }
        }
    }
    return;
}

void sorted_trie::trie_node::find_within(std::deque<int>& collected, const std::vector<char>& current, int offset, int limit, size_t iter) {
    if (!scores) {
        scores = new std::vector<int>(current.size()+1);
    } else {
        scores->resize(current.size() + 1);
    }

    // Initializing the first level of the levenshtein array here. We multiply the limit by MULT to handle fractions.
    {
        int score=0;
        for (auto& s : *scores) {
            s=score;
            score+=EDIT;
        }
        limit*=MULT;
    }

    if (children) {
        for (size_t i=0; i<NBASES; ++i) {
            auto& child=(*children)[i];
            if (!child.dead_end()) {
                child.find_within(collected, current, offset, BASES[i], *scores, limit, iter);
            }
        }
    }

    return;
}

/****************************************
* Exposed class for users of the trie. *
****************************************/

sorted_trie::sorted_trie(const size_t nseq, const char** seqs, const int* lens) : counter(0) {
    current.reserve(50); // probably enough.
    for (size_t i=0; i<nseq; ++i) {
        const char * cur_str=seqs[i];
        const size_t cur_len=lens[i];
        current.resize(cur_len);
        std::copy(cur_str, cur_str+cur_len, current.begin());
        toplevel.insert(i, current);
    }
    current.clear();
    return;
}

const std::deque<int>& sorted_trie::find(const char* cur_str, const int cur_len, int limit) {
    // Figuring out the common prefix from the last 'current'.
    size_t idx=0;
    while (idx < current.size() && idx < cur_len && cur_str[idx]==current[idx]) {
        ++idx;
    }
    const size_t common=idx;

    // If it's fully common, this implies that this string is the same as the previous string,
    // so we use the previous results and skip to avoid redundant work.
    if (common==cur_len && cur_len==current.size() && counter) {
        return collected;
    }

    // Replacing the remaining elements in 'current'.
    current.resize(cur_len);
    while (idx < cur_len) {
        current[idx]=cur_str[idx];
        ++idx;
    }

//    for (size_t i=0; i<cur_len; ++i) {
//        Rprintf("%c", current[i]);
//    }
//    Rprintf("\n");
//    Rprintf("Common is %i\n", common);

    // Searching through the trie.
    toplevel.find_within(collected, current, common, limit, counter);

    ++counter;
    return collected;
}

void sorted_trie::dump() {
    toplevel.dump();
    return;
}

void sorted_trie::order(size_t n, const char ** seqs, const int* lens, int* order) {
	std::iota(order, order+n, 0);
	std::sort(order, order+n, [&] (int left, int right) -> bool {
		const char * L=seqs[left], *R=seqs[right];
		const int LN=lens[left], RN=lens[right], limit=std::min(LN, RN);
		for (int i=0; i<limit; ++i) {
			if (L[i]!=R[i]) { return L[i] < R[i]; }
		}
		return LN < RN;
	});
    return;
}


