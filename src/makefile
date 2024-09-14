# Define variables
CXX = g++
CXXFLAGS = -pthread -lz
CPP_DIR = cpp
TARGETS = statistic_test_combination kmerpos
SRC_STAT = $(CPP_DIR)/statistic_test_combination.cpp $(CPP_DIR)/statistic_test_main.cpp
SRC_KMER = $(CPP_DIR)/kmerpos.cpp $(CPP_DIR)/main.cpp

# Default target
all: $(TARGETS)

# Generate statistic_test_combination
statistic_test_combination: $(SRC_STAT)
	@echo "Compiling statistic_test_combination..."
	$(CXX) -o $@ $^ $(CXXFLAGS)

# Generate kmerpos
kmerpos: $(SRC_KMER)
	@echo "Compiling kmerpos..."
	$(CXX) -o $@ $^ $(CXXFLAGS)

# Clean target
clean:
	@echo "Cleaning up..."
	rm -f $(TARGETS)

.PHONY: all clean
