
CC=g++
PARAMS=-O3 -Wall
BUILD_DIR=build

all: $(BUILD_DIR)/all-hole-gradient $(BUILD_DIR)/all-hole-hessian $(BUILD_DIR)/contraction-test

$(BUILD_DIR)/all-hole-gradient: all-hole-gradient.cpp contraction/contraction.h
	$(CC) $(PARAMS) -fopenmp -o $(BUILD_DIR)/all-hole-gradient all-hole-gradient.cpp

$(BUILD_DIR)/all-hole-hessian: all-hole-hessian.cpp contraction/contraction.h
	$(CC) $(PARAMS) -fopenmp -o $(BUILD_DIR)/all-hole-hessian all-hole-hessian.cpp

$(BUILD_DIR)/contraction-test: contraction-test.cpp contraction/contraction.h test/test-data.h
	$(CC) $(PARAMS) -fopenmp -o $(BUILD_DIR)/contraction-test contraction-test.cpp -lgtest -lgtest_main

clean: $(BUILD_DIR)/all-hole-gradient $(BUILD_DIR)/all-hole-hessian $(BUILD_DIR)/contraction-test
	rm $(BUILD_DIR)/all-hole-gradient
	rm $(BUILD_DIR)/all-hole-hessian
	rm $(BUILD_DIR)/contraction-test
