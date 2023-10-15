# High Performance Contraction of Brickwall Quantum Circuits for Hamiltonian Simulation

This software allows contracting brickwall tensor networks efficiently, obtaining their gradient and Hessian. It uses gate-wise application and caching to reduce memory usage to <2GB on 16 qubits. Performance is improved by a factor 10-30 compared to a simple Kronecker product based matrix contraction. For more details see my master thesis.

## Building

Make sure you have `make` and a modern version of the GCC C++ compiler `g++` installed. Then run:

```
make
```

to build the two sample programs and the unit tests.

## Running

To execute the unit tests after building, run:

```
./build/contraction-test
```

To use the computation as part of your own application, take a look at `all-hole-gradient.cpp` and `all-hole-hessian.cpp` as an example on how to use the contraction code as a blackbox for computing gradient and Hessian of your brickwall circuits.

## Configuration
The contraction process is configured using a `CircuitGateConfiguration`. It can be generated from a list of pairs specifying the qubits on which the gates act using the `generateLayerConfigurationFromGates()` function from `contraction/contraction_utility.h`.

You also might want to take a look at changing `PARALLEL_THREAD_COUNT` in `contraction/contraction.h` to match your hardware.