CXX ?= g++
CXXFLAGS ?= -O2 -Wall -std=c++20 -pthread -Isrc

SRC := src
BIN := bin

.PHONY: all build check clean

all: $(BIN)/test_xortools $(BIN)/full_anneal $(BIN)/random_regular

check: $(BIN)/test_xortools
	./$(BIN)/test_xortools

$(BIN):
	mkdir -p $(BIN)

$(BIN)/test_xortools: $(BIN)/test_xortools.o $(BIN)/bitmatrix.o $(BIN)/utils.o | $(BIN)
	$(CXX) $(CXXFLAGS) $(BIN)/test_xortools.o $(BIN)/bitmatrix.o $(BIN)/utils.o -o $@

$(BIN)/full_anneal: $(BIN)/full_anneal.o $(BIN)/bitmatrix.o $(BIN)/utils.o $(BIN)/xoropt.o | $(BIN)
	$(CXX) $(CXXFLAGS) $(BIN)/full_anneal.o $(BIN)/bitmatrix.o $(BIN)/utils.o $(BIN)/xoropt.o -o $@

$(BIN)/random_regular: $(BIN)/random_regular.o $(BIN)/bitmatrix.o $(BIN)/matrix_ensembles.o $(BIN)/utils.o | $(BIN)
	$(CXX) $(CXXFLAGS) $(BIN)/random_regular.o $(BIN)/bitmatrix.o $(BIN)/matrix_ensembles.o $(BIN)/utils.o -o $@

$(BIN)/test_xortools.o: $(SRC)/test_xortools.cpp $(SRC)/bitmatrix.hpp $(SRC)/utils.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN)/bitmatrix.o: $(SRC)/bitmatrix.cpp $(SRC)/bitmatrix.hpp $(SRC)/utils.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN)/full_anneal.o: $(SRC)/full_anneal.cpp $(SRC)/bitmatrix.hpp $(SRC)/utils.hpp $(SRC)/xoropt.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN)/random_regular.o: $(SRC)/random_regular.cpp $(SRC)/bitmatrix.hpp $(SRC)/matrix_ensembles.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN)/xoropt.o: $(SRC)/xoropt.cpp $(SRC)/bitmatrix.hpp $(SRC)/utils.hpp $(SRC)/xoropt.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN)/matrix_ensembles.o: $(SRC)/matrix_ensembles.cpp $(SRC)/bitmatrix.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(BIN)/utils.o: $(SRC)/utils.cpp $(SRC)/utils.hpp | $(BIN)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -rf $(BIN)
	rm -f data/densetest.txt data/sparsetest.txt data/vec_densetest.txt data/vec_sparsetest.txt
