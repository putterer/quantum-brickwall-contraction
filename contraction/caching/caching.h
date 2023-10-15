#ifndef CONTRACTION_CACHING_H
#define CONTRACTION_CACHING_H

#include <complex>

using namespace std;

// The indexing allows picking the correct state vectors for layer 0 (from the right) by accessing left[0] and right[0]
// -> the state vectors required for computing gate[2] can be found in left[2] and right[2]

// 256 MB of statevectors
template <int L, int Q, int S> struct ContractionCache {
    static const int layers = L;
    static const int qubits = Q;
    static const int statevecSize = S;

    complex<double> right[L][Q/2][S]; // indexing starts at the right -> layers already applied
    complex<double> left[L][Q/2][S]; // indexing starts at the right, 1 offset -> layers still left to be applied
};


#endif