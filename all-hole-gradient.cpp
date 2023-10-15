#include <stdio.h>
#include <cstdlib>

#define QUBITS 12
#define STATEVEC_SIZE (1 << QUBITS)

#define LAYERS 8

#define CACHED_CONTRACTION

#include "contraction/contraction.h"

using namespace std;


int main() {

    CircuitGateConfiguration config = generate1DBrickwallNonOffset<LAYERS, QUBITS>();

    // 6 KiB on stack
    complex<double> gates[LAYERS][4][4];
    complex<double> gradients[LAYERS][4][4];

    // initialize gates
    for(int layer = 0;layer < LAYERS;layer++) {
        for(int row = 0;row < 4;row++) {
            for(int column = 0;column < 4;column++) {
                gates[layer][row][column] = randComplex();
                gradients[layer][row][column] = complex(0.0, 0.0);
            }
        }
    }

    // Calculate gradient for each layer = for each hole
#ifndef CACHED_CONTRACTION
    computeAllFullGradients<LAYERS, QUBITS, STATEVEC_SIZE, false>(gates, gradients, config);
#endif

#ifdef CACHED_CONTRACTION
    computeAllFullGradients<LAYERS, QUBITS, STATEVEC_SIZE, true>(gates, gradients, config);
#endif

    cout << gradients[0][1][1] << endl;
}