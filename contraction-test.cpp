#include "gtest/gtest.h"


#define QUBITS 12
#define STATEVEC_SIZE (1 << QUBITS)

#include "contraction/contraction.h"
#include "contraction/contraction_utility.h"


#include "test/test-data.h"


using namespace std;

inline void expectEqualsStatevec(complex<double> *b1, complex<double> *b2, int statevecSize, double orderOfMagnitudeThreshold = -10) {
    for(int i = 0;i < statevecSize;i++) {
        if(b1[i].real() != b2[i].real()) {
            EXPECT_LE(log10f(abs(b1[i].real() - b2[i].real()) / abs(b1[i].real())), orderOfMagnitudeThreshold) << "Statevecs differ at index " << i << " " << b1[i] << " vs " << b2[i];
        }
        if(b1[i].imag() != b2[i].imag()) {
            EXPECT_LE(log10f(abs(b1[i].imag() - b2[i].imag()) / abs(b1[i].imag())), orderOfMagnitudeThreshold) << "Statevecs differ at index " << i << " " << b1[i] << " vs " << b2[i];
        }
    }
}

inline void expectEqualsMatrix(complex<double> m1[][4], complex<double> m2[][4], int rows, int cols, double orderOfMagnitudeThreshold = -10) {
    for(int row = 0;row < rows;row++) {
        expectEqualsStatevec(m1[row], m2[row], cols, orderOfMagnitudeThreshold);
    }
}

template <int layers, int qubits, int statevecSize> complex<double> contractFullNetwork(complex<double> gates[layers][4][4], int oneHotIndex, CircuitGateConfiguration<layers, qubits> &config) {
    complex<double> *buffer1 = new complex<double>[STATEVEC_SIZE];
    complex<double> *buffer2 = new complex<double>[STATEVEC_SIZE];

    complex<double> *inputBuffer = buffer1;
    complex<double> *outputBuffer = buffer2;

    for(int i = 0;i < STATEVEC_SIZE;i++) {
        inputBuffer[i] = complex(0.0, 0.0);
    }

    inputBuffer[oneHotIndex] = complex(1.0, 0.0);

    for(int layer = 0;layer < layers;layer++) {
        applyGateAtAllSites<qubits, statevecSize>(&inputBuffer, &outputBuffer, gates[layer], config.layerConfig[layer]);
    }

    complex<double> res = inputBuffer[oneHotIndex];

    delete [] buffer1;
    delete [] buffer2;

    // vector contraction without hole = dot product
    // 1-hot vec -> just pick from other
    return res;
}


TEST(SingleSiteApplication, preservesIdentity) {

    complex<double> *inputBuffer = new complex<double>[STATEVEC_SIZE];
    complex<double> *outputBuffer = new complex<double>[STATEVEC_SIZE];

    for(int i = 0;i < STATEVEC_SIZE;i++) {
        inputBuffer[i] = randComplex();
    }

    // identity gate
    complex<double> gate[4][4] = {
        1.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i,
        0.0 + 0.0i, 1.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i,
        0.0 + 0.0i, 0.0 + 0.0i, 1.0 + 0.0i, 0.0 + 0.0i,
        0.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i, 1.0 + 0.0i
    };

    applySingleGate<QUBITS, STATEVEC_SIZE>(2, 3, inputBuffer, outputBuffer, gate);
    
    expectEqualsStatevec(inputBuffer, outputBuffer, STATEVEC_SIZE);

    delete [] inputBuffer;
    delete [] outputBuffer;
}

TEST(AllSiteApplication, preservesIdentity) {

    LayerGateConfiguration layerConfig = generateWrappingLayer<QUBITS>(0);

    complex<double> *buffer1 = new complex<double>[STATEVEC_SIZE];
    complex<double> *buffer2 = new complex<double>[STATEVEC_SIZE];

    complex<double> *inputBuffer = buffer1;
    complex<double> *outputBuffer = buffer2;

    for(int i = 0;i < STATEVEC_SIZE;i++) {
        inputBuffer[i] = randComplex();
    }

    // identity gate
    complex<double> gate[4][4] = {
        1.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i,
        0.0 + 0.0i, 1.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i,
        0.0 + 0.0i, 0.0 + 0.0i, 1.0 + 0.0i, 0.0 + 0.0i,
        0.0 + 0.0i, 0.0 + 0.0i, 0.0 + 0.0i, 1.0 + 0.0i
    };

    applyGateAtAllSites<QUBITS, STATEVEC_SIZE>(&inputBuffer, &outputBuffer, gate, layerConfig);
    
    expectEqualsStatevec(inputBuffer, outputBuffer, STATEVEC_SIZE);

    delete [] buffer1;
    delete [] buffer2;
}



template <int qubits> inline void checkTestCase(TestData<qubits> &testData) {
    LayerGateConfiguration layerConfig = generateWrappingLayer<qubits>(0);

    complex<double> *output = new complex<double>[testData.statevecSize];

    complex<double> *inputBuffer = testData.input; // this is going to overwrite the testData.input array! Don't reuse!
    complex<double> *outputBuffer = output;

    complex<double> gate[4][4];
    memcpy(gate, testData.gate, sizeof(gate));

    for(int i = 0;i < testData.iterations;i++) {
        applyGateAtAllSites<testData.qubits, testData.statevecSize>(&inputBuffer, &outputBuffer, gate, layerConfig);
    }
    
    printf("\n");
    printf("Actual:  \t%.25lf,\t%.25lf\n", inputBuffer[0].real(), inputBuffer[0].imag());
    printf("Expected:\t%.25lf,\t%.25lf\n", testData.output[0].real(), testData.output[0].imag());
    printf("Error:   \t%.25lf,\t%.25lf\n", abs(inputBuffer[0].real() - testData.output[0].real()), abs(inputBuffer[0].imag() - testData.output[0].imag()));
    double numericalError = abs(inputBuffer[0].real() - testData.output[0].real());
    cout << "Numerical error: " << (numericalError / inputBuffer[0].real()) << "    Order of Magnitude: " << log10f(abs(numericalError / inputBuffer[0].real())) << endl;
    printf("\n");

    expectEqualsStatevec(inputBuffer, testData.output, testData.statevecSize);

    delete [] output;
}

TEST(AllSiteApplication, randomGateOneHotInput) {
    checkTestCase(TEST_RANDOM_GATE_ONE_HOT);
}

TEST(AllSiteApplication, randomGateRandomInput) {
    checkTestCase(TEST_RANDOM_GATE_RANDOM_VEC);
}

TEST(FourLayerApplication, randomGateRandomInput) {
    checkTestCase(TEST_RANDOM_GATE_RANDOM_VEC_4);
}

TEST(LargeEightLayerApplication, randomGateRandomInput) {
    checkTestCase(TEST_LARGE_RANDOM_GATE_RANDOM_VEC_8);
}

TEST(HoleContraction, symmetricExample) {
    complex<double> vec1[8] = {
        0.0 + 0.0i, 1.0 + 0.0i, 2.0 + 0.0i, 3.0 + 0.0i, 4.0 + 0.0i, 5.0 + 0.0i, 6.0 + 0.0i, 7.0 + 0.0i
    };
    complex<double> vec2[8] = {
        0.0 + 0.0i, 1.0 + 0.0i, 2.0 + 0.0i, 3.0 + 0.0i, 4.0 + 0.0i, 5.0 + 0.0i, 6.0 + 0.0i, 7.0 + 0.0i
    };

    complex<double> res[4][4];

    contractHoleStatevecs<3, 8>(vec1, vec2, 0, 1, res);

    for(int i = 0;i < 4;i++) {
        for(int j = 0;j < 4;j++) {
            cout << res[i][j];
        }
        cout << endl;
    }

    complex<double> expected[4][4] = {
        16.0 + 0.0i, 20.0 + 0.0i, 24.0 + 0.0i, 28.0 + 0.0i,
        20.0 + 0.0i, 26.0 + 0.0i, 32.0 + 0.0i, 38.0 + 0.0i,
        24.0 + 0.0i, 32.0 + 0.0i, 40.0 + 0.0i, 48.0 + 0.0i,
        28.0 + 0.0i, 38.0 + 0.0i, 48.0 + 0.0i, 58.0 + 0.0i
    };

    expectEqualsMatrix(expected, res, 4, 4);
}

TEST(HoleContraction, nonSymmetricExample) {
    complex<double> vec1[8] = {
        0.0 + 0.0i, 1.0 + 0.0i, 2.0 + 0.0i, 3.0 + 0.0i, 4.0 + 0.0i, 5.0 + 0.0i, 6.0 + 0.0i, 7.0 + 0.0i
    };
    complex<double> vec2[8] = {
        1.0 + 0.0i, 1.0 + 0.0i, 2.0 + 0.0i, 3.0 + 0.0i, 4.0 + 0.0i, 5.0 + 0.0i, 6.0 + 0.0i, 7.0 + 0.0i
    };

    complex<double> res[4][4];

    contractHoleStatevecs<3, 8>(vec1, vec2, 0, 1, res);

    for(int i = 0;i < 4;i++) {
        for(int j = 0;j < 4;j++) {
            cout << res[i][j];
        }
        cout << endl;
    }

    complex<double> expected[4][4] = {
        16.0 + 0.0i, 20.0 + 0.0i, 24.0 + 0.0i, 28.0 + 0.0i,
        21.0 + 0.0i, 26.0 + 0.0i, 32.0 + 0.0i, 38.0 + 0.0i,
        26.0 + 0.0i, 32.0 + 0.0i, 40.0 + 0.0i, 48.0 + 0.0i,
        31.0 + 0.0i, 38.0 + 0.0i, 48.0 + 0.0i, 58.0 + 0.0i
    };

    expectEqualsMatrix(expected, res, 4, 4);
}

TEST(Gradient, finiteDiffApprox) {
    const int layers = 4;
    const int qubits = 8;
    const int statevecSize = 256;

    CircuitGateConfiguration config = generate1DBrickwallNonOffset<layers, qubits>();

    complex<double> gates[layers][4][4];
    complex<double> gradients[layers][4][4];
    complex<double> gradientsEggspected[layers][4][4];

    for(int layer = 0;layer < layers;layer++) {
        for(int row = 0;row < 4;row++) {
            for(int column = 0;column < 4;column++) {
                gates[layer][row][column] = randComplex();
                gradients[layer][row][column] = complex(0.0, 0.0);
                gradientsEggspected[layer][row][column] = complex(0.0, 0.0);
            }
        }
    }

    // UUT
    computeAllFullGradients<layers, qubits, statevecSize>(gates, gradients, config);

    // approximate expected result using a finite difference approximation
    complex<double> delta = 1.0E-8 + 0.0i;

    for(int oneHotIndex = 0;oneHotIndex < statevecSize;oneHotIndex++) {
        for(int gate = 0;gate < layers;gate++) {
            for(int i = 0;i < 4;i++) {
                for(int j = 0;j < 4;j++) {
                    gates[gate][i][j] += delta;
                    complex<double> higher = contractFullNetwork<layers, qubits, statevecSize>(gates, oneHotIndex, config);
                    gates[gate][i][j] -= 2.0 * delta;
                    complex<double> lower = contractFullNetwork<layers, qubits, statevecSize>(gates, oneHotIndex, config);
                    gates[gate][i][j] += delta;

                    complex<double> paramGrad = (higher - lower) / delta / 2.0;

                    gradientsEggspected[gate][i][j] += paramGrad;
                }
            }
        }
    }

    cout << endl << "Exact hole computation:" << endl;

    for(int i = 0;i < 4;i++) {
        for(int j = 0;j < 4;j++) {
            cout << gradients[0][i][j];
        }
        cout << endl;
    }

    cout << endl << "Finite difference approximation:" << endl;

    for(int i = 0;i < 4;i++) {
        for(int j = 0;j < 4;j++) {
            cout << gradientsEggspected[0][i][j];
        }
        cout << endl;
    }

    cout << endl;

    for(int layer = 0;layer < layers;layer++) {
        expectEqualsMatrix(gradientsEggspected[layer], gradients[layer], 4, 4, -3.0);
    }
}

TEST(Hessian, finiteDiffApprox) {
    const int layers = 4;
    const int qubits = 6;
    const int statevecSize = 64;

    CircuitGateConfiguration config = generate1DBrickwallNonOffset<layers, qubits>();

    complex<double> gates[layers][4][4];
    complex<double> (*hessians)[layers][layers][4][4][4][4] = new complex<double>[1][layers][layers][4][4][4][4];
    complex<double> (*hessiansEggspected)[layers][layers][4][4][4][4] = new complex<double>[1][layers][layers][4][4][4][4];

    for(int layer = 0;layer < layers;layer++) {
        for(int row = 0;row < 4;row++) {
            for(int column = 0;column < 4;column++) {
                gates[layer][row][column] = randComplex();
            }
        }
    }

    for(int layer1 = 0;layer1 < layers;layer1++) {
        for(int layer2 = 0;layer2 < layers;layer2++) {
            for(int row1 = 0;row1 < 4;row1++) {
                for(int column1 = 0;column1 < 4;column1++) {
                    for(int row2 = 0;row2 < 4;row2++) {
                        for(int column2 = 0;column2 < 4;column2++) {
                            hessians[0][layer1][layer2][row1][column1][row2][column2] = complex(0.0, 0.0);
                            hessiansEggspected[0][layer1][layer2][row1][column1][row2][column2] = complex(0.0, 0.0);
                        }
                    }
                }
            }
        }
    }

    // UUT
    computeAllFullHessians<layers, qubits, statevecSize>(gates, hessians[0], config);

    // approximate expected result using a finite difference approximation
    complex<double> delta = 1.0E-4 + 0.0i;

    for(int oneHotIndex = 0;oneHotIndex < statevecSize;oneHotIndex++) {
        
        cout << oneHotIndex << ", " << flush;

        for(int layer1 = 0;layer1 < layers;layer1++) {
            for(int layer2 = 0;layer2 < layers;layer2++) {
                for(int row1 = 0;row1 < 4;row1++) {
                    for(int column1 = 0;column1 < 4;column1++) {
                        for(int row2 = 0;row2 < 4;row2++) {
                            for(int column2 = 0;column2 < 4;column2++) {
                                complex<double> hess = 0.0;

                                gates[layer1][row1][column1] += delta;
                                gates[layer2][row2][column2] += delta; // 1+, 2+
                                hess += contractFullNetwork<layers, qubits, statevecSize>(gates, oneHotIndex, config);
                                gates[layer2][row2][column2] -= 2.0 * delta; // 1+, 2-
                                hess -= contractFullNetwork<layers, qubits, statevecSize>(gates, oneHotIndex, config);
                                gates[layer1][row1][column1] -= 2.0 * delta; // 1-, 2-
                                hess += contractFullNetwork<layers, qubits, statevecSize>(gates, oneHotIndex, config);
                                gates[layer2][row2][column2] += 2.0 * delta; // 1-, 2+
                                hess -= contractFullNetwork<layers, qubits, statevecSize>(gates, oneHotIndex, config);
                                gates[layer1][row1][column1] += delta; // 1
                                gates[layer2][row2][column2] -= delta; // 2
                                
                                hess /= 4.0 * delta * delta;

                                hessiansEggspected[0][layer1][layer2][row1][column1][row2][column2] += hess;
                            }
                        }
                    }
                }
            }
        }
    }

    for(int layer1 = 0;layer1 < layers;layer1++) {
        for(int layer2 = 0;layer2 < layers;layer2++) {
            for(int row1 = 0;row1 < 4;row1++) {
                for(int column1 = 0;column1 < 4;column1++) {
                    expectEqualsMatrix(hessiansEggspected[0][layer1][layer2][row1][column1], hessians[0][layer1][layer2][row1][column1], 4, 4, -3.0);
                }
            }
        }
    }

    cout << endl << "Exact hole computation:" << endl;

    for(int i = 0;i < 4;i++) {
        for(int j = 0;j < 4;j++) {
            cout << hessians[0][2][2][3][3][i][j];  // THE FIRST INDEX IS NOT!!! AN INDEX! (it's a pointer, not an array)
        }
        cout << endl;
    }

    cout << endl << "Finite difference approximation:" << endl;

    for(int i = 0;i < 4;i++) {
        for(int j = 0;j < 4;j++) {
            cout << hessiansEggspected[0][2][2][3][3][i][j];  // THE FIRST INDEX IS NOT!!! AN INDEX!
        }
        cout << endl;
    }

    cout << endl;
}

TEST(Splitting, oneHotMatrices) {
    const int qubits = 12;
    const int statevecSize = (1 << qubits);

    complex<double> *inputStatevec = new complex<double>[statevecSize];

    for(int i = 0;i < statevecSize;i++) {
        inputStatevec[i] = randComplex();
    }

    StatevecBufferPair<statevecSize> parallelStatevecs[4][4];
    StatevecBufferPair<statevecSize> expectedParallelStatevecs[4][4];

    // UUT
    computeHoleSplit<qubits, statevecSize>(inputStatevec, parallelStatevecs, 2, 3);

    for(int row = 0;row < 4;row++) {
        for(int column = 0;column < 4;column++) {
            complex<double> oneHotMat[4][4];
            for(int i = 0;i < 4;i++) {
                for(int j = 0;j < 4;j++) {
                    oneHotMat[i][j] = 0;
                }
            }

            oneHotMat[row][column] = 1;

            applySingleGate<qubits, statevecSize>(2, 3, inputStatevec, expectedParallelStatevecs[row][column].inputBuffer, oneHotMat);
        }
    }

    for(int row = 0;row < 4;row++) {
        for(int column = 0;column < 4;column++) {
            expectEqualsStatevec(expectedParallelStatevecs[row][column].inputBuffer, parallelStatevecs[row][column].inputBuffer, statevecSize);
        }
    }

    for(int i = 0;i < 10;i++) {
        cout << parallelStatevecs[1][2].inputBuffer[i] << " ";
    }
    cout << endl;
    for(int i = 0;i < 10;i++) {
        cout << expectedParallelStatevecs[1][2].inputBuffer[i] << " ";
    }
    cout << endl;

    delete inputStatevec;
}
