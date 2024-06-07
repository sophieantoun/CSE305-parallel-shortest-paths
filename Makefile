CXX = g++
CXXFLAGS = -pthread -std=c++11 -Wall

# Executable name
EXECUTABLE = output_program

# Source files
SOURCES = main.cpp DeltaSteppingSequential.cpp delta_stepping_parallel.cpp Dijkstra.cpp 

# Object files
OBJECTS = $(SOURCES:.cpp=.o)

# Rule to link the executable
$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS)

# Rule to compile object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean workspace
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

# Declare phony targets
.PHONY: clean
