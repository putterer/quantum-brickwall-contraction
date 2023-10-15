#ifndef CONTRACTION_H
#define CONTRACTION_H

#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <string.h>
#include <omp.h>

#include <complex>

#include "gate_application.h"
#include "contraction_utility.h"
#include "caching/caching.h"

using namespace std;

#define PARALLEL_THREAD_COUNT 12


template <int layers, int qubits, int statevecSize>
void generateContractionCache(complex<double> gates[layers][4][4], complex<double> gatesInverted[layers][4][4], int oneHotIndex, ContractionCache<layers, qubits, statevecSize> *targetCache);


template <int qubits, int statevecSize>
void contractHoleStatevecs(complex<double> *leftGateOutStatevec, complex<double> *rightGateInStatevec, int q1, int q2, complex<double> (&gradientTarget)[4][4]) {
    // gradient[i][j] = sum(others) { rightIn[others | j] * leftOut[others | i] }

    // this is the same as a single gate application
    // instead of multiplying with a gate, we are multiplying with a differently indexed 2nd statevector

    int smallerQubits = std::min(q1, q2);
    int centerQubits = abs(q2 - q1) - 1;
    int gateQubits = 2;
    int biggerQubits = qubits - gateQubits - smallerQubits - centerQubits;

    for(int biggerQubitsState = 0;biggerQubitsState < (1 << biggerQubits);biggerQubitsState++) { // slower running, bigger endianess
        for(int centerQubitsState = 0;centerQubitsState < (1 << centerQubits);centerQubitsState++) {
            for(int smallerQubitsState = 0;smallerQubitsState < (1 << smallerQubits);smallerQubitsState++) {  // the faster running, smaller endianess qubits
                int otherQubitsState = smallerQubitsState | (centerQubitsState << (smallerQubits + 1)) | (biggerQubitsState << (smallerQubits + centerQubits + gateQubits));

                for(int q2StateLeft = 0;q2StateLeft < (1 << 1);q2StateLeft++) {
                    for(int q1StateLeft = 0;q1StateLeft < (1 << 1);q1StateLeft++) {
                        int leftEntryIndex = otherQubitsState | (q2StateLeft << q2) | (q1StateLeft << q1);
                        int holeQubitsLeftGateOut = q1StateLeft | (q2StateLeft << 1);

                        for(int q2StateRight = 0;q2StateRight < (1 << 1);q2StateRight++) {
                            for(int q1StateRight = 0;q1StateRight < (1 << 1);q1StateRight++) {
                                int rightEntryIndex = otherQubitsState | (q2StateRight << q2) | (q1StateRight << q1);
                                int holeQubitsRightGateIn = q1StateRight | (q2StateRight << 1);

                                gradientTarget[holeQubitsLeftGateOut][holeQubitsRightGateIn] += leftGateOutStatevec[leftEntryIndex] * rightGateInStatevec[rightEntryIndex];
                            }
                        }
                    }
                }
            }
        }
    }
}

// Gates / layers are sorted from right to left
template <int layers, int qubits, int statevecSize, bool cachingEnabled>
void computeHoleGradient(complex<double> gates[layers][4][4], complex<double> gatesInverted[layers][4][4], int holeLayer, int holeQubit, int oneHotIndex, complex<double> (&gradientTarget)[4][4], ContractionCache<layers, qubits, statevecSize> *cache, CircuitGateConfiguration<layers, qubits> &config) {
    
    if(cachingEnabled) {
        int cacheSiteIndex = config.layerConfig[holeLayer].sitesByLowerQubit[holeQubit];

        complex<double> *rightInputBuffer = cache->right[holeLayer][cacheSiteIndex];
        complex<double> *leftInputBuffer = cache->left[holeLayer][cacheSiteIndex];

        int secondQubit = config.layerConfig[holeLayer].mapping[holeQubit];
        contractHoleStatevecs<qubits, statevecSize>(leftInputBuffer, rightInputBuffer, holeQubit, secondQubit, gradientTarget);

    } else {
        complex<double> *leftBuffer1 = new complex<double>[statevecSize];
        complex<double> *leftBuffer2 = new complex<double>[statevecSize];
        complex<double> *rightBuffer1 = new complex<double>[statevecSize];
        complex<double> *rightBuffer2 = new complex<double>[statevecSize];

        complex<double> *leftInputBuffer = leftBuffer1;
        complex<double> *rightInputBuffer = rightBuffer1;
        complex<double> *leftOutputBuffer = leftBuffer2;
        complex<double> *rightOutputBuffer = rightBuffer2;

        // initialize buffers
        for(int i = 0;i < statevecSize;i++) {
            leftInputBuffer[i] = complex(0.0, 0.0);
            rightInputBuffer[i] = complex(0.0, 0.0);
        }

        leftInputBuffer[oneHotIndex] = complex(1.0, 0.0);
        rightInputBuffer[oneHotIndex] = complex(1.0, 0.0);

        // this applies a 1-hot vector on both sides
        // if you want to apply a matrix operator, you need to load its row here and in the cache generation

        // Contract both sides of the circuit
        // Right
        for(int layer = 0;layer < holeLayer;layer++) {
            applyGateAtAllSites<qubits, statevecSize>(&rightInputBuffer, &rightOutputBuffer, gates[layer], config.layerConfig[layer]);
        }

        // Left
        for(int layer = layers - 1;layer > holeLayer;layer--) {
            // this order would be incorrect for gate wise caching! (which is off here)
            applyGateAtAllSites<qubits, statevecSize>(&leftInputBuffer, &leftOutputBuffer, gatesInverted[layer], config.layerConfig[layer]);
        }

        // Hole layer
        applyGateAtAllNonHoleSites<qubits, statevecSize>(&rightInputBuffer, &rightOutputBuffer, gates[holeLayer], holeQubit, config.layerConfig[holeLayer]);

        int secondQubit = config.layerConfig[holeLayer].mapping[holeQubit];
        contractHoleStatevecs<qubits, statevecSize>(leftInputBuffer, rightInputBuffer, holeQubit, secondQubit, gradientTarget);

        delete [] leftBuffer1;
        delete [] leftBuffer2;
        delete [] rightBuffer1;
        delete [] rightBuffer2;
    }

}


template <int layers, int qubits, int statevecSize, bool cachingEnabled>
void computeGateGradient(complex<double> gates[layers][4][4], complex<double> gatesInverted[layers][4][4], int gateLayer, int oneHotIndex, complex<double> (&gradientTarget)[4][4], ContractionCache<layers, qubits, statevecSize> *cache, CircuitGateConfiguration<layers, qubits> &config) {
    for(int gateHoleIndex : config.layerConfig[gateLayer].lowerQubitsBySite) {
        computeHoleGradient<layers, qubits, statevecSize, cachingEnabled>(gates, gatesInverted, gateLayer, gateHoleIndex, oneHotIndex, gradientTarget, cache, config);
    }
}


template <int layers, int qubits, int statevecSize, bool cachingEnabled>
void computeAllGradients(complex<double> gates[layers][4][4], complex<double> gatesInverted[layers][4][4], int oneHotIndex, complex<double> (&gradients)[layers][4][4], ContractionCache<layers, qubits, statevecSize> *cache, CircuitGateConfiguration<layers, qubits> &config) {
    for(int layer = 0;layer < layers;layer++) {
        // for each layer: compute gradient of each gate (hole) -> sum them up
        computeGateGradient<layers, qubits, statevecSize, cachingEnabled>(gates, gatesInverted, layer, oneHotIndex, gradients[layer], cache, config);
    }
}

template <int layers, int qubits, int statevecSize, bool cachingEnabled>
void computeAllGradients(complex<double> gates[layers][4][4], int oneHotIndex, complex<double> (&gradients)[layers][4][4], ContractionCache<layers, qubits, statevecSize> *cache, CircuitGateConfiguration<layers, qubits> &config) {
    complex<double> gatesInverted[layers][4][4];
    invertGates<layers>(gates, gatesInverted);
    
    computeAllGradients<layers, qubits, statevecSize, cachingEnabled>(gates, gatesInverted, oneHotIndex, gradients, cache, config);
}

template <int layers, int qubits, int statevecSize, bool cachingEnabled = true>
void computeAllFullGradients(complex<double> gates[layers][4][4], complex<double> (&gradients)[layers][4][4], CircuitGateConfiguration<layers, qubits> &config) {
    complex<double> gatesInverted[layers][4][4];
    invertGates<layers>(gates, gatesInverted);

    int threads = PARALLEL_THREAD_COUNT;
    // allocates 2 KiB * 12 = 24 KiB on the stack
    complex<double> parallelGradients[threads][layers][4][4];
    
    // initialize gradients
    for(int thread = 0;thread < threads;thread++) {
        for(int layer = 0;layer < layers;layer++) {
            for(int row = 0;row < 4;row++) {
                for(int column = 0;column < 4;column++) {
                    parallelGradients[thread][layer][row][column] = complex(0.0, 0.0);
                }
            }
        }
    }


    #pragma omp parallel for num_threads(threads)
    for(int oneHotIndex = 0;oneHotIndex < statevecSize;oneHotIndex++) {
        int threadId = omp_get_thread_num();

        ContractionCache<layers, qubits, statevecSize> *cache = nullptr;
        
        if(cachingEnabled) {
            cache = new ContractionCache<layers, qubits, statevecSize>();
            generateContractionCache<layers, qubits, statevecSize>(gates, gatesInverted, oneHotIndex, cache, config);
        }

        computeAllGradients<layers, qubits, statevecSize, cachingEnabled>(gates, gatesInverted, oneHotIndex, parallelGradients[threadId], cache, config);

        if(cachingEnabled) {
            delete cache;
        }
    }

    // reduce parallel results
    for(int thread = 0;thread < threads;thread++) {
        for(int layer = 0;layer < layers;layer++) {
            for(int row = 0;row < 4;row++) {
                for(int column = 0;column < 4;column++) {
                    gradients[layer][row][column] += parallelGradients[thread][layer][row][column];
                }
            }
        }
    }
}



// -----------------------------
// Hessian
// -----------------------------

template <int statevecSize>
struct StatevecBufferPair
{
    complex<double> *buffer1 = new complex<double>[statevecSize];
    complex<double> *buffer2 = new complex<double>[statevecSize];
    complex<double> *inputBuffer = buffer1;
    complex<double> *outputBuffer = buffer2;

    ~StatevecBufferPair() {
        delete [] buffer1;
        delete [] buffer2;
    }
};

template <int qubits, int statevecSize>
void computeHoleSplit(complex<double> *inputStatevec, StatevecBufferPair<statevecSize> (&outputStatevecs)[4][4], int holeSiteQ1, int holeSiteQ2) {

    // I hope you like bit shift logic :)

    int holeQ1IndexMask = 0b1 << holeSiteQ1;
    int holeQ2IndexMask = 0b1 << holeSiteQ2;
    int holeIndexMask = holeQ1IndexMask | holeQ2IndexMask;
    int otherIndexMask = 0b1111111111111111 ^ holeIndexMask; // this only works for up to 16 qubits

    for(int inputIndex = 0;inputIndex < statevecSize;inputIndex++) {
        int column = (inputIndex & holeQ1IndexMask) >> holeSiteQ1 | ((inputIndex & holeQ2IndexMask) >> holeSiteQ2) << 1;
        int other = inputIndex & otherIndexMask;

        for(int row = 0;row < 4;row++) {
            int outputIndex = other | ((row & 0b01) << holeSiteQ1) | (((row & 0b10) >> 1) << holeSiteQ2);
            outputStatevecs[row][column].inputBuffer[outputIndex] = inputStatevec[inputIndex]; // only one column, one other affects it -> equality, no summation necessary
        }
    }
}


template <int layers, int qubits, int statevecSize>
void computeHessianWithoutCache(complex<double> gates[layers][4][4], complex<double> gatesInverted[layers][4][4], int holeLayer1, int holeSite1, int holeLayer2, int holeSite2, int oneHotIndex, complex<double> (&hessian)[4][4][4][4], CircuitGateConfiguration<layers, qubits> &config) {
    complex<double> *leftBuffer1 = new complex<double>[statevecSize];
    complex<double> *leftBuffer2 = new complex<double>[statevecSize];
    complex<double> *rightBuffer1 = new complex<double>[statevecSize];
    complex<double> *rightBuffer2 = new complex<double>[statevecSize];

    complex<double> *leftInputBuffer = leftBuffer1;
    complex<double> *leftOutputBuffer = leftBuffer2;
    complex<double> *rightInputBuffer = rightBuffer1;
    complex<double> *rightOutputBuffer = rightBuffer2;

    // run right before hole1
    for(int layer = 0;layer < holeLayer1;layer++) {
        applyGateAtAllSites<qubits, statevecSize>(&rightInputBuffer, &rightOutputBuffer, gates[layer], config.layerConfig[layer]);
    }
    applyGateAtAllSitesBefore<qubits, statevecSize>(&rightInputBuffer, &rightOutputBuffer, gates[holeLayer1], holeSite1, config.layerConfig[holeLayer1]);

    // run left after hole2
    for(int layer = layers - 1;layer > holeLayer2;layer--) {
        // this order would be incorrect for gate wise caching! (which is off here)
        applyGateAtAllSites<qubits, statevecSize>(&leftInputBuffer, &leftOutputBuffer, gatesInverted[layer], config.layerConfig[layer]);
    }
    applyGateAtAllSitesAfter<qubits, statevecSize>(&leftInputBuffer, &leftOutputBuffer, gatesInverted[holeLayer2], holeSite2, config.layerConfig[holeLayer2]);


    // split at hole1

    // allocate parallel statevecs, 16 MB in total, the buffers themselves are allocated on the heap
    StatevecBufferPair<statevecSize> parallelStatevecs[4][4];

    int splitSecondQubit = config.layerConfig[holeLayer1].mapping[holeSite1];
    computeHoleSplit<qubits, statevecSize>(rightInputBuffer, parallelStatevecs, holeSite1, splitSecondQubit);


    // run center, hole1 -> hole2, for every statevec
    for(int row = 0;row < 4;row++) {
        for(int column = 0;column < 4;column++) {
            applyGateAtAllSitesAfter<qubits, statevecSize>(&(parallelStatevecs[row][column].inputBuffer), &(parallelStatevecs[row][column].outputBuffer), gates[holeLayer1], holeSite1, config.layerConfig[holeLayer1]);

            for(int layer = holeLayer1 + 1;layer < holeLayer2;layer++) {
                applyGateAtAllSites<qubits, statevecSize>(&(parallelStatevecs[row][column].inputBuffer), &(parallelStatevecs[row][column].outputBuffer), gates[layer], config.layerConfig[layer]);
            }

            applyGateAtAllSitesBefore<qubits, statevecSize>(&(parallelStatevecs[row][column].inputBuffer), &(parallelStatevecs[row][column].outputBuffer), gates[holeLayer2], holeSite2, config.layerConfig[holeLayer2]);
        }
    }

    
    // contract all at hole2
    for(int row = 0;row < 4;row++) {
        for(int column = 0;column < 4;column++) {
            int secondQubit = config.layerConfig[holeLayer2].mapping[holeSite2];
            contractHoleStatevecs<qubits, statevecSize>(leftInputBuffer, parallelStatevecs[row][column].inputBuffer, holeSite2, secondQubit, hessian[row][column]);
        }
    }


    delete [] leftBuffer1;
    delete [] leftBuffer2;
    delete [] rightBuffer1;
    delete [] rightBuffer2;
}

template <int layers, int qubits, int statevecSize>
void computeAllHessiansWithoutCache(complex<double> gates[layers][4][4], complex<double> gatesInverted[layers][4][4], int oneHotIndex, complex<double> (&hessians)[layers][layers][4][4][4][4], CircuitGateConfiguration<layers, qubits> &config) {
    
    // hole1 is always before hole2

    for(int holeLayer1 = 0;holeLayer1 < layers;holeLayer1++) {
        for(int holeSite1 : config.layerConfig[holeLayer1].lowerQubitsBySite) {
            if(holeLayer1 == layers - 1 && holeSite1 == config.layerConfig[holeLayer1].lastSite) {
                continue; // if this is the last hole, we can't find a second hole
            }

            for(int holeLayer2 = holeLayer1;holeLayer2 < layers;holeLayer2++) {
                for(int holeSite2 : config.layerConfig[holeLayer2].lowerQubitsBySite) {
                    if(holeLayer1 == holeLayer2 && holeSite2 <= holeSite1) {
                        continue;
                    }

                    computeHessianWithoutCache<layers, qubits, statevecSize>(gates, gatesInverted, holeLayer1, holeSite1, holeLayer2, holeSite2, oneHotIndex, hessians[holeLayer1][holeLayer2], config);
                }
            }
        }
    }
}




template <int layers, int qubits, int statevecSize>
void computeAllHessiansForFirstHoleCached(complex<double> gates[layers][4][4], complex<double> gatesInverted[layers][4][4], int holeLayer1, int holeSite1, int oneHotIndex, complex<double> (&hessians)[layers][4][4][4][4], ContractionCache<layers, qubits, statevecSize> *cache, CircuitGateConfiguration<layers, qubits> &config) {
    
    int cacheSiteIndexRight = config.layerConfig[holeLayer1].sitesByLowerQubit[holeSite1];
    complex<double> *rightCached = cache->right[holeLayer1][cacheSiteIndexRight];

    // split at hole1

    // allocate parallel statevecs, 16 MB in total, the buffers themselves are allocated on the heap
    StatevecBufferPair<statevecSize> parallelStatevecs[4][4];

    int splitSecondQubit = config.layerConfig[holeLayer1].mapping[holeSite1];
    computeHoleSplit<qubits, statevecSize>(rightCached, parallelStatevecs, holeSite1, splitSecondQubit);


    // run right till end
    int currentLayer = holeLayer1; // indicates next site

    bool applying = false;

    for(;currentLayer < layers;currentLayer++) {
        for(int currentSite : config.layerConfig[currentLayer].lowerQubitsBySite) {
            if(currentLayer == holeLayer1 && currentSite == holeSite1) {
                applying = true;
                continue;
            }
            if(!applying) {
                continue;
            }

            // contract site
            int cacheSiteIndexCurrent = config.layerConfig[currentLayer].sitesByLowerQubit[currentSite];
            complex<double> *leftCached = cache->left[currentLayer][cacheSiteIndexCurrent];

            for(int row = 0;row < 4;row++) {
                for(int column = 0;column < 4;column++) {
                    int secondQubit = config.layerConfig[currentLayer].mapping[currentSite];
                    contractHoleStatevecs<qubits, statevecSize>(leftCached, parallelStatevecs[row][column].inputBuffer, currentSite, secondQubit, hessians[currentLayer][row][column]);
                }
            }
            
            
            // run site, all 16 vecs
            for(int row = 0;row < 4;row++) {
                for(int column = 0;column < 4;column++) {
                    applySingleGateWithSwap<qubits, statevecSize>(currentSite, &(parallelStatevecs[row][column].inputBuffer), &(parallelStatevecs[row][column].outputBuffer), gates[currentLayer], config.layerConfig[currentLayer]);
                }
            }

            // this unnecessarily runs for the last index, can that be removed in reasonable complexity?

        }
    }

}

template <int layers, int qubits, int statevecSize>
void computeAllHessiansCached(complex<double> gates[layers][4][4], complex<double> gatesInverted[layers][4][4], int oneHotIndex, complex<double> (&hessians)[layers][layers][4][4][4][4], CircuitGateConfiguration<layers, qubits> &config) {
    
    // run right, run left, full, -> cache
    ContractionCache<layers, qubits, statevecSize> *cache = new ContractionCache<layers, qubits, statevecSize>();
    generateContractionCache<layers, qubits, statevecSize>(gates, gatesInverted, oneHotIndex, cache, config);
    
    for(int holeLayer1 = 0;holeLayer1 < layers;holeLayer1++) {
        for(int holeSite1 : config.layerConfig[holeLayer1].lowerQubitsBySite) {
            computeAllHessiansForFirstHoleCached<layers, qubits, statevecSize>(gates, gatesInverted, holeLayer1, holeSite1, oneHotIndex, hessians[holeLayer1], cache, config);
        }
    }

    delete cache;
}



template <int layers, int qubits, int statevecSize, bool cachingEnabled = true>
void computeAllFullHessians(complex<double> gates[layers][4][4], complex<double> (&hessians)[layers][layers][4][4][4][4], CircuitGateConfiguration<layers, qubits> &config) {

    complex<double> gatesInverted[layers][4][4];
    invertGates<layers>(gates, gatesInverted);

    int threads = PARALLEL_THREAD_COUNT;
    
    // allocates 3.2 MB
    complex<double> (*parallelHessians)[layers][layers][4][4][4][4] = new complex<double>[threads][layers][layers][4][4][4][4];

    // initialize hessians
    for(int thread = 0;thread < threads;thread++) {
        for(int layer1 = 0;layer1 < layers;layer1++) {
            for(int layer2 = 0;layer2 < layers;layer2++) {
                for(int row1 = 0;row1 < 4;row1++) {
                    for(int column1 = 0;column1 < 4;column1++) {
                        for(int row2 = 0;row2 < 4;row2++) {
                            for(int column2 = 0;column2 < 4;column2++) {
                                parallelHessians[thread][layer1][layer2][row1][column1][row2][column2] = complex(0.0, 0.0);
                            }
                        }
                    }
                }
            }
        }
    }


    #pragma omp parallel for num_threads(threads)
    for(int oneHotIndex = 0;oneHotIndex < statevecSize;oneHotIndex++) {
        int threadId = omp_get_thread_num();

        if(cachingEnabled) {
            computeAllHessiansCached<layers, qubits, statevecSize>(gates, gatesInverted, oneHotIndex, parallelHessians[threadId], config);
        } else {
            computeAllHessiansWithoutCache<layers, qubits, statevecSize>(gates, gatesInverted, oneHotIndex, parallelHessians[threadId], config);
        }
        

    }


    // reduce parallel results
    for(int thread = 0;thread < threads;thread++) {
        for(int layer1 = 0;layer1 < layers;layer1++) {
            for(int layer2 = 0;layer2 < layers;layer2++) {
                for(int row1 = 0;row1 < 4;row1++) {
                    for(int column1 = 0;column1 < 4;column1++) {
                        for(int row2 = 0;row2 < 4;row2++) {
                            for(int column2 = 0;column2 < 4;column2++) {
                                hessians[layer1][layer2][row1][column1][row2][column2] += parallelHessians[thread][layer1][layer2][row1][column1][row2][column2];
                            }
                        }
                    }
                }
            }
        }
    }

    // mirror hessian tensor
    for(int layer1 = 0;layer1 < layers;layer1++) {
        for(int layer2 = 0;layer2 < layer1;layer2++) {
            for(int row1 = 0;row1 < 4;row1++) {
                for(int column1 = 0;column1 < 4;column1++) {
                    for(int row2 = 0;row2 < 4;row2++) {
                        for(int column2 = 0;column2 < 4;column2++) {
                            hessians[layer1][layer2][row1][column1][row2][column2] = hessians[layer2][layer1][row2][column2][row1][column1];
                        }
                    }
                }
            }
        }
    }

    // duplicate same layer
    for(int layer = 0;layer < layers;layer++) {
        for(int row1 = 0;row1 < 4;row1++) {
            for(int column1 = 0;column1 < 4;column1++) {
                for(int row2 = 0;row2 < 4;row2++) {
                    for(int column2 = 0;column2 < 4;column2++) {
                        hessians[layer][layer][row1][column1][row2][column2] *= 2.0;
                    }
                }
            }
        }
    }

    delete [] parallelHessians;
}





// ------------------------------
// Caching
// ------------------------------

template <int layers, int qubits, int statevecSize>
void generateContractionCache(complex<double> gates[layers][4][4], complex<double> gatesInverted[layers][4][4], int oneHotIndex, ContractionCache<layers, qubits, statevecSize> *targetCache, CircuitGateConfiguration<layers, qubits> &config) {

    // initialize buffers, 1 MiB each
    complex<double> *buffer1 = new complex<double>[statevecSize];
    complex<double> *buffer2 = new complex<double>[statevecSize];

    complex<double> *inputBuffer = buffer1;
    complex<double> *outputBuffer = buffer2;

    for(int i = 0;i < statevecSize;i++) {
        inputBuffer[i] = complex(0.0, 0.0);
    }

    inputBuffer[oneHotIndex] = complex(1.0, 0.0);


    // layer indicates what gate will be applied next
    // cache layer is the layer that would be applied next

    // run circuit right -> left
    for(int layer = 0;layer < layers;layer++) {
        // apply layer
        for(int siteIndex = 0;siteIndex < qubits / 2;siteIndex++) {
            int applicationSite = config.layerConfig[layer].lowerQubitsBySite[siteIndex];

            for(int i = 0;i < statevecSize;i++) {
                targetCache->right[layer][siteIndex][i] = inputBuffer[i];
            }

            int secondQubit = config.layerConfig[layer].mapping[applicationSite];
            applySingleGate<qubits, statevecSize>(applicationSite, secondQubit, inputBuffer, outputBuffer, gates[layer]);
        
            // swap
            complex<double> *tmp = inputBuffer;
            inputBuffer = outputBuffer;
            outputBuffer = tmp;
        }
    }


    // reinitialize buffer
    for(int i = 0;i < statevecSize;i++) {
        inputBuffer[i] = complex(0.0, 0.0);
    }

    inputBuffer[oneHotIndex] = complex(1.0, 0.0);
    
    // this applies a 1-hot vector on both sides
    // if you want to apply a matrix operator, you need to load its row here and in the hole contraction

    // run circuit left -> right
    for(int layer = layers - 1;layer >= 0;layer--) {
        // apply layer
        for(int siteIndex = (qubits / 2) - 1;siteIndex >= 0;siteIndex--) {
            int applicationSite = config.layerConfig[layer].lowerQubitsBySite[siteIndex];

            for(int i = 0;i < statevecSize;i++) {
                targetCache->left[layer][siteIndex][i] = inputBuffer[i];
            }

            int secondQubit = config.layerConfig[layer].mapping[applicationSite];
            applySingleGate<qubits, statevecSize>(applicationSite, secondQubit, inputBuffer, outputBuffer, gatesInverted[layer]);
        
            // swap
            complex<double> *tmp = inputBuffer;
            inputBuffer = outputBuffer;
            outputBuffer = tmp;
        }
    }

    // free buffers
    delete [] buffer1;
    delete [] buffer2;
}


#endif