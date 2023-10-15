#include <stdio.h>
#include <cstdlib>

#define QUBITS 8
#define STATEVEC_SIZE (1 << QUBITS)

#define LAYERS 8

#define CACHED_CONTRACTION

#include "contraction/contraction.h"

using namespace std;


int main() {

    CircuitGateConfiguration config = generate1DBrickwallNonOffset<LAYERS, QUBITS>();

    // 6 KiB on stack
    complex<double> gates[LAYERS][4][4];
    complex<double> (*hessians)[LAYERS][LAYERS][4][4][4][4] = new complex<double>[1][LAYERS][LAYERS][4][4][4][4];

    // initialize gates
    for(int layer = 0;layer < LAYERS;layer++) {
        for(int row = 0;row < 4;row++) {
            for(int column = 0;column < 4;column++) {
                gates[layer][row][column] = randComplex();
            }
        }
    }

    // initialize hessians
    for(int layer1 = 0;layer1 < LAYERS;layer1++) {
        for(int layer2 = 0;layer2 < LAYERS;layer2++) {
            for(int row1 = 0;row1 < 4;row1++) {
                for(int column1 = 0;column1 < 4;column1++) {
                    for(int row2 = 0;row2 < 4;row2++) {
                        for(int column2 = 0;column2 < 4;column2++) {
                            hessians[0][layer1][layer2][row1][column1][row2][column2] = complex(0.0, 0.0);
                        }
                    }
                }
            }
        }
    }

    // Calculate gradient for each layer = for each hole
#ifndef CACHED_CONTRACTION
    computeAllFullHessians<LAYERS, QUBITS, STATEVEC_SIZE, false>(gates, hessians[0], config);
#endif

#ifdef CACHED_CONTRACTION
    computeAllFullHessians<LAYERS, QUBITS, STATEVEC_SIZE, true>(gates, hessians[0], config);
#endif

    cout << hessians[0][3][1][1][2][3][3] << endl;

    delete [] hessians;
}
