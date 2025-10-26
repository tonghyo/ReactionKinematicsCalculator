# Makefile for FusionReaction simulation
# Requires ROOT installation

# ROOT configuration
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++11 -Wall -O2 $(ROOTCFLAGS)
LIBS = $(ROOTLIBS)

# Source files
SOURCES = FusionReaction_Setup.cpp FusionReaction_MassHist.cpp FusionReaction_Kinematics.cpp FusionReaction_Analysis.cpp
HEADERS = FusionReaction.h
MAIN = fusion_reaction.C

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Executable
TARGET = fusion_reaction

# Default target
all: $(TARGET)

# Build executable
$(TARGET): $(OBJECTS) $(MAIN)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LIBS)

# Build object files
%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build files
clean:
	rm -f $(OBJECTS) $(TARGET) *.root *.png *.pdf

# Run the simulation
run: $(TARGET)
	./$(TARGET)

# Help
help:
	@echo "Available targets:"
	@echo "  all     - Build the executable"
	@echo "  clean   - Remove build files"
	@echo "  run     - Build and run the simulation"
	@echo "  help    - Show this help message"

.PHONY: all clean run help
