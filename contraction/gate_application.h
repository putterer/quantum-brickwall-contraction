#ifndef GATE_APPLCATION_H
#define GATE_APPLCATION_H

#include <stdio.h>
#include <complex>
#include <algorithm>
#include <stdlib.h>

#include "contraction_utility.h"

using namespace std;

template <int qubits, int statevecSize> void applySingleGate(int q1, int q2, complex<double> inputStatevec[statevecSize], complex<double> outputStatevec[statevecSize], complex<double> gate[4][4]) {
    int smallerQubits = std::min(q1, q2);
    int centerQubits = abs(q2 - q1) - 1;
    int gateQubits = 2;
    int biggerQubits = qubits - gateQubits - smallerQubits - centerQubits;

    // for all other qubits
    for(int biggerQubitsState = 0;biggerQubitsState < (1 << biggerQubits);biggerQubitsState++) { // slower running, bigger endianess
        for(int centerQubitsState = 0;centerQubitsState < (1 << centerQubits);centerQubitsState++) {
            for(int smallerQubitsState = 0;smallerQubitsState < (1 << smallerQubits);smallerQubitsState++) {  // the faster running, smaller endianess qubits
                int otherQubitsState = smallerQubitsState | (centerQubitsState << (smallerQubits + 1)) | (biggerQubitsState << (smallerQubits + centerQubits + gateQubits));

                for(int q2StateOutput = 0;q2StateOutput < (1 << 1);q2StateOutput++) {
                    for(int q1StateOutput = 0;q1StateOutput < (1 << 1);q1StateOutput++) {
                        int outputEntryIndex = otherQubitsState | (q2StateOutput << q2) | (q1StateOutput << q1);
                        int gateQubitsStateOutput = q1StateOutput | (q2StateOutput << 1);

                        complex<double> outputEntry = 0.0;
    
                        for(int q2StateInput = 0;q2StateInput < (1 << 1);q2StateInput++) {
                            for(int q1StateInput = 0;q1StateInput < (1 << 1);q1StateInput++) {
                                int inputEntryIndex = otherQubitsState | (q2StateInput << q2) | (q1StateInput << q1);
                                int gateQubitsStateInput = q1StateInput | (q2StateInput << 1);

                                complex<double> gateFactor = gate[gateQubitsStateOutput][gateQubitsStateInput]; // how to cache?

                                complex<double> inputEntry = inputStatevec[inputEntryIndex];

                                outputEntry += inputEntry * gateFactor;
                            }
                        }
    
                        outputStatevec[outputEntryIndex] = outputEntry;
                    }
                }
            }
        }
    }
}


template <int qubits, int statevecSize> void applySingleGateWithSwap(int lowerQubitToApplyTo, complex<double> **inputStatevec, complex<double> **outputStatevec, complex<double> gate[4][4], LayerGateConfiguration<qubits> &layerConfig) {
    int secondQubit = layerConfig.mapping[lowerQubitToApplyTo];
    applySingleGate<qubits, statevecSize>(lowerQubitToApplyTo, secondQubit, *inputStatevec, *outputStatevec, gate);

    // swap
    complex<double> *tmp = *inputStatevec;
    *inputStatevec = *outputStatevec;
    *outputStatevec = tmp;
}

// takes the input and output buffer as a pass-by-pointer array, so they can be swapped in this function
template <int qubits, int statevecSize> void applyGateAtAllSites(complex<double> **inputStatevec, complex<double> **outputStatevec, complex<double> gate[4][4], LayerGateConfiguration<qubits> &layerConfig) {
    for(int applicationSite : layerConfig.lowerQubitsBySite) {
        int secondQubit = layerConfig.mapping[applicationSite];
        applySingleGate<qubits, statevecSize>(applicationSite, secondQubit, *inputStatevec, *outputStatevec, gate);
        
        // swap
        complex<double> *tmp = *inputStatevec;
        *inputStatevec = *outputStatevec;
        *outputStatevec = tmp;
    }
}

// Applies the gate at all not excluded sites
// this could have been solved using a default template argument, but would lead to less clear code and dependence on C++17, so the function is duplicated
template <int qubits, int statevecSize> void applyGateAtAllNonHoleSites(complex<double> **inputStatevec, complex<double> **outputStatevec, complex<double> gate[4][4], int excludedApplicationSite, LayerGateConfiguration<qubits> &layerConfig) {
    for(int applicationSite : layerConfig.lowerQubitsBySite) {
        if(applicationSite == excludedApplicationSite) {
            continue;
        }

        int secondQubit = layerConfig.mapping[applicationSite];
        applySingleGate<qubits, statevecSize>(applicationSite, secondQubit, *inputStatevec, *outputStatevec, gate);
        
        // swap
        complex<double> *tmp = *inputStatevec;
        *inputStatevec = *outputStatevec;
        *outputStatevec = tmp;
    }
}

// Applies the gate at all sites before param
template <int qubits, int statevecSize> void applyGateAtAllSitesBefore(complex<double> **inputStatevec, complex<double> **outputStatevec, complex<double> gate[4][4], int firstExcludedSite, LayerGateConfiguration<qubits> &layerConfig) {
    for(int applicationSite : layerConfig.lowerQubitsBySite) {
        if(applicationSite == firstExcludedSite) {
            break;
        }

        int secondQubit = layerConfig.mapping[applicationSite];
        applySingleGate<qubits, statevecSize>(applicationSite, secondQubit, *inputStatevec, *outputStatevec, gate);
        
        // swap
        complex<double> *tmp = *inputStatevec;
        *inputStatevec = *outputStatevec;
        *outputStatevec = tmp;
    }
}

// Applies the gate at all sites after param
template <int qubits, int statevecSize> void applyGateAtAllSitesAfter(complex<double> **inputStatevec, complex<double> **outputStatevec, complex<double> gate[4][4], int firstExcludedSite, LayerGateConfiguration<qubits> &layerConfig) {
    bool applying = false;
    for(int applicationSite : layerConfig.lowerQubitsBySite) {
        if(applicationSite == firstExcludedSite) {
            applying = true;
            continue;
        }
        if(!applying) {
            continue;
        }

        int secondQubit = layerConfig.mapping[applicationSite];
        applySingleGate<qubits, statevecSize>(applicationSite, secondQubit, *inputStatevec, *outputStatevec, gate);
        
        // swap
        complex<double> *tmp = *inputStatevec;
        *inputStatevec = *outputStatevec;
        *outputStatevec = tmp;
    }
}


#endif