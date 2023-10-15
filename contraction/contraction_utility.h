#ifndef CONTRACTION_UTILITY_H
#define CONTRACTION_UTILITY_H

#include <cstdlib>
#include "complex.h"

using namespace std;

#define MOD(x, m) ((x % m + m) % m)


template <int qubits> struct LayerGateConfiguration
{
    int mapping[qubits];  // specifies a mapping from each qubit to the other qubit on its gate
    int lowerQubitsBySite[qubits / 2];  // specifies a list of qubits to iterate over to find all gates
    int sitesByLowerQubit[qubits];  // reverse lookupTable
    int lastSite;  // lowerQubitsBySite[-1] for fast lookup
};

// describes the orientation of gates in the circuit
template <int layers, int qubits> struct CircuitGateConfiguration
{
    struct LayerGateConfiguration<qubits> layerConfig[layers];
};

// returns a layer configuration corresponding to a list of qubit pairs
template <int qubits> struct LayerGateConfiguration<qubits> generateLayerConfigurationFromGates(int gateQubits[qubits / 2][2]) {
    struct LayerGateConfiguration<qubits> res;

    for(int gate = 0;gate < qubits / 2;gate++) {
        res.mapping[gateQubits[gate][0]] = gateQubits[gate][1];
        res.mapping[gateQubits[gate][1]] = gateQubits[gate][0];
        res.lowerQubitsBySite[gate] = gateQubits[gate][0];
        res.sitesByLowerQubit[gateQubits[gate][0]] = gate;
    }

    res.lastSite = gateQubits[qubits / 2 - 1][0];

    return res;
}

// generates a 1D layer with an offset, that wraps if necessary
template <int qubits> struct LayerGateConfiguration<qubits> generateWrappingLayer(int offset) {
    struct LayerGateConfiguration<qubits> res;
    for(int q = 0;q < qubits;q++) {
        if(q % 2 == 0) {
            res.mapping[MOD(q + offset, qubits)] = MOD(q + 1 + offset, qubits);
            res.lowerQubitsBySite[q/2] = MOD(q + offset, qubits);
            res.sitesByLowerQubit[MOD(q + offset, qubits)] = q/2;
            res.lastSite = res.lowerQubitsBySite[q/2];
        } else {
            res.mapping[MOD(q + offset, qubits)] = MOD(q - 1 + offset, qubits);
        }
    }
    return res;
}

// generates a normal 1D brickwall network from wrapping layers
template <int layers, int qubits> struct CircuitGateConfiguration<layers, qubits> generate1DBrickwall() {
    struct CircuitGateConfiguration<layers, qubits> res;
    for(int l = 0;l < layers;l++) {
        res.layerConfig[l] = generateWrappingLayer<qubits>(l % 2);
    }
    return res;
}

// generates a 1D brickwall network without any wrap around or offsets, for testing purposes
template <int layers, int qubits> struct CircuitGateConfiguration<layers, qubits> generate1DBrickwallNonOffset() {
    struct CircuitGateConfiguration<layers, qubits> res;
    for(int l = 0;l < layers;l++) {
        res.layerConfig[l] = generateWrappingLayer<qubits>(0);
    }
    return res;
}



void invertGate(complex<double> gate[4][4], complex<double> (&gateInvertedTarget)[4][4]) {
    for(int i = 0;i < 4;i++) {
        for(int j = 0;j < 4;j++) {
            gateInvertedTarget[j][i] = gate[i][j];
        }
    }
}

template <int layers> void invertGates(complex<double> gates[layers][4][4], complex<double> (&gatesInvertedTarget)[layers][4][4]) {
    for(int layer = 0;layer < layers;layer++) {
        invertGate(gates[layer], gatesInvertedTarget[layer]);
    }
}



double randIntveral(double min, double max) {
    double uniform = ((double)rand()) / ((double)RAND_MAX);
    return uniform * (max - min) + min;
}

complex<double> randComplex() {
    complex<double> res(randIntveral(-1.0, 1.0), randIntveral(-1.0, 1.0));
    return res;
}

template <int bits> void printBinary(int x) {
    for(int i = 0;i < bits;i++) {
        if(x & 1) {
            printf("1");
        } else {
            printf("0");
        }
        x >>= 1;
    }
}

#endif