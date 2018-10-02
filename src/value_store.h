#ifndef VALUE_STORE_H
#define VALUE_STORE_H
#include "sarlacc.h"

/* This defines a class that stores multiple variable length arrays in a single container.
 * This avoids the need to reallocate memory for separate deques/vectors for different arrays
 * or for loops where a similar amount of memory needs to be used across iterations.
 *
 * Arrays can be of variable length but are assumed to not change in size after creation.
 */ 

template<typename T, class V>
class value_store {
public:
    value_store(size_t n=100, size_t n2=1000) : storage_values(n2), starts(n), lengths(n), n_stored(0), nvals_stored(0) {}

    template<class IT>
    void add(IT start, IT end, size_t index) {
        add_internal(start, end, index);
        if (index >= n_stored) {
            if (index > n_stored) {
                std::fill(lengths.begin() + n_stored, lengths.begin() + index, 0); // Clearing intervening values to ensure validity.
            }
            n_stored=index + 1;
        }
        return;
    }

    template<class IT>
    void add(IT start, IT end) {
        add_internal(start, end, n_stored);
        ++n_stored;
        return;
    }
    
    void clear() {
        n_stored=0;
        nvals_stored=0;
        return;
    }

    typename V::const_iterator get_start(size_t i) const {
        return storage_values.begin() + starts[i];
    }

    typename V::iterator get_start_unsafe(size_t i) {
        return storage_values.begin() + starts[i];
    }

    const size_t get_len(size_t i) const {
        return lengths[i];
    }

    const size_t size () const {
        return n_stored;
    }
private:
    V storage_values;
    std::deque<size_t> starts, lengths;
    size_t n_stored, nvals_stored;

    template<class IT>
    void add_internal(IT start, IT end, size_t index) {
        if (index >= starts.size()) {
            starts.resize(index+1);
            lengths.resize(index+1);
        } 

        const auto len=end - start;
        starts[index]=nvals_stored;
        lengths[index]=len;

        const auto required=nvals_stored + len;
        if (required > storage_values.size()) {
            storage_values.resize(required);
        }

        std::copy(start, end, storage_values.begin() + nvals_stored);
        nvals_stored=required;
        return;
    }
};


#endif
