#ifndef SORTED_TRIE_H
#define SORTED_TRIE_H
#include <vector>
#include <deque>
#include <algorithm>

class sorted_trie {
public:
    sorted_trie(const size_t, const char**, const int*);
    const std::deque<int>& find(const char*, const int, int);
    void dump();

    static void order(const size_t, const char**, const int*, int*);
private:
    std::vector<char> current;
    std::deque<int> collected;
    size_t counter;

    struct trie_node {
        std::deque<int>* indices;
        std::vector<trie_node>* children;
        std::vector<int>* scores;
        size_t history;
        
        trie_node() {}
        ~trie_node() {}

        bool dead_end () const {}
        void insert (const int, const std::vector<char>&, size_t);
        void insert(int, const std::vector<char>&);
        
        void dump(int, char);
        void dump();

        void find_within(std::deque<int>&, const std::vector<char>&, size_t, char, std::vector<int>&, int, size_t);
        void find_within(std::deque<int>&, const std::vector<char>&, int, int, size_t);
    };

    trie_node toplevel;
};

#endif
