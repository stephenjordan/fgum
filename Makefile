CXX ?= g++
CPPFLAGS ?= -Isrc
CXXFLAGS ?= -O2 -Wall -std=c++20 -pthread
LDFLAGS ?=
LDLIBS ?=

SRC := src
BIN := bin

PROGRAMS := \
	$(BIN)/test_xortools \
	$(BIN)/full_anneal \
	$(BIN)/random_regular \
	$(BIN)/density

.PHONY: all build check clean

all build: $(PROGRAMS)

check: $(BIN)/test_xortools $(BIN)/density
	./$(BIN)/test_xortools
	./$(BIN)/density --k 2 --D 3 --bins 101 --iterations 2 --precision 0.1 --max-delta 0.1 >/dev/null

$(BIN):
	mkdir -p $(BIN)

$(BIN)/test_xortools: $(BIN)/test_xortools.o $(BIN)/bitmatrix.o $(BIN)/utils.o | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BIN)/full_anneal: $(BIN)/full_anneal.o $(BIN)/bitmatrix.o $(BIN)/utils.o $(BIN)/xoropt.o | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BIN)/random_regular: $(BIN)/random_regular.o $(BIN)/bitmatrix.o $(BIN)/utils.o | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BIN)/density: $(BIN)/density.o $(BIN)/density_evolution.o $(BIN)/distribution.o $(BIN)/utils.o | $(BIN)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

$(BIN)/%.o: $(SRC)/%.cpp | $(BIN)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(wildcard $(BIN)/*.d)

clean:
	rm -rf $(BIN)
	rm -f data/densetest.txt data/sparsetest.txt data/vec_densetest.txt data/vec_sparsetest.txt
