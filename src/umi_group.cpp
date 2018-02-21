#include "sarlacc.h"

/**********************************************************************
 * A functor to calculate the levenshtein distance between DNAStrings.
 **********************************************************************/
 
struct lev_dist {
private:
    const XStringSet_holder * ptr;
public:
    lev_dist(const XStringSet_holder* p) : ptr(p) {};

    int operator()(const int left, const int right) const {
        auto lData=get_elt_from_XStringSet_holder(ptr, left);
        auto rData=get_elt_from_XStringSet_holder(ptr, right);
        const char * sL=lData.ptr, * sR=rData.ptr;
        const size_t lenL = lData.length, lenR = rData.length;

        // Modified from http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance
        std::vector<std::vector<unsigned int> > d(lenL + 1, std::vector<unsigned int>(lenR + 1));
        
        d[0][0] = 0;
        for(unsigned int i = 1; i <= lenL; ++i) {
            d[i][0] = i;
        }
        for(unsigned int i = 1; i <= lenR; ++i) {
            d[0][i] = i;
        }
        
        for (unsigned int i = 1; i <= lenL; ++i) {
            for(unsigned int j = 1; j <= lenR; ++j) {
                // note that std::min({arg1, arg2, arg3}) works only in C++11,
                d[i][j] = std::min({ d[i - 1][j] + 1, d[i][j - 1] + 1, d[i - 1][j - 1] + (DNAdecode(sL[i - 1]) == DNAdecode(sR[j - 1]) ? 0 : 1) });
            }
        }

        return d[lenL][lenR]; 
    }
};

/*
 * BK-tree implementation in C++
 * Copyright (C) 2012 Eiichi Sato
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* This file was obtained from https://gist.github.com/eiiches/2016232, 2018/02/21.
 * I have modified it to use the Biostrings access methods.
 */

#include <map>
#include <cmath>
#include <vector>

namespace qq {

typedef int MetricType;
typedef int KeyType;

/**********************************************************************
 * A class for the nodes of the BK tree.
 **********************************************************************/

namespace detail {

class tree_node {
private:
	KeyType value;
	std::map<MetricType, tree_node*> * children;
    const lev_dist * distptr;

public:
	tree_node(const KeyType &key, const lev_dist& d) : value(key), children(NULL), distptr(&d) { }

	~tree_node() {
		if (children) {
			for (auto iter = children->begin(); iter != children->end(); ++iter) {
				delete iter->second;
            }
			delete children;
		}
	}

	bool insert(tree_node *node) {
		if (!node) {
			return false;
        }
	
        MetricType distance = (*distptr)(node->value, this->value);
		if (distance == 0) {
			return false; /* value already exists */
        }

		if (!children) {
			children = new std::map<MetricType, tree_node*>();
        }

		auto iterator = children->find(distance);
		if (iterator == children->end()) {
			children->insert(std::make_pair(distance, node));
			return true;
		}

		return iterator->second->insert(node);
	}

protected:
	bool has_children() const {
		return this->children && this->children->size();
	}

	void _find_within(std::vector<std::pair<KeyType, MetricType>> &result, const KeyType &key, MetricType d) const {
		MetricType n = (*distptr)(key, this->value);
        if (n <= d) {
			result.push_back(std::make_pair(this->value, n));
        }

		if (!this->has_children()) {
			return;
        }

        // Triangle inequality magic applies here.
		for (auto iter = children->begin(); iter != children->end(); ++iter) {
			MetricType distance = iter->first;
			if (n - d <= distance && distance <= n + d) {
				iter->second->_find_within(result, key, d);
            }
		}
	}

public:
	std::vector<std::pair<KeyType, MetricType>> find_within(const KeyType &key, MetricType d) const {
		std::vector<std::pair<KeyType, MetricType>> result;
		_find_within(result, key, d);
		return result;
	}

	void dump_tree(int depth = 0) {
		for (int i = 0; i < depth; ++i) {
			Rprintf("    ");
        }
		Rprintf("%i\n", this->value);
		if (this->has_children()) {
            for (auto iter = children->begin(); iter != children->end(); ++iter) {
                iter->second->dump_tree(depth + 1);
            }
        }
	}
};

} /* namespace detail */

class bktree {
private:
    detail::tree_node *m_top;
	size_t m_n_nodes;
    lev_dist dist_fun;

public:
	bktree(const XStringSet_holder* p) : m_top(NULL), m_n_nodes(0), dist_fun(p) { }

public:
	void insert(const KeyType &key) {
        detail::tree_node *node = new detail::tree_node(key, dist_fun);
		if (!m_top) {
			m_top = node;
			m_n_nodes = 1;
			return;
		}
		if (m_top->insert(node)) {
			++m_n_nodes;
        }
	};

public:
	std::vector<std::pair<KeyType, MetricType>> find_within(KeyType key, MetricType d) const {
		return m_top->find_within(key, d);
	}

	void dump_tree() {
		m_top->dump_tree();
	}

public:
	size_t size() const {
		return m_n_nodes;
	}
};

} /* namespace qq */


SEXP umi_group (SEXP umi, SEXP threshold) {
    BEGIN_RCPP
        
    // Checking inputs.       
    auto seq=hold_XStringSet(umi); 
    const size_t nseq=get_length_from_XStringSet_holder(&seq);
    int limit=check_integer_scalar(threshold, "distance threshold");

    // Creating the BK-tree.
    qq::bktree dic(&seq);
    for (size_t i=0; i<nseq; ++i) {
        Rprintf("Currently at: %i\n", i);
        dic.insert(i);
        auto out=dic.find_within(i, limit);
        Rprintf("%i\n", out.size());
    }
   
    return Rcpp::IntegerVector::create(1);
    END_RCPP
}
