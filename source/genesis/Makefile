# Genesis Protein Evolution Project Makefile

# Compiler and flags
CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2
DEBUG_FLAGS = -std=c++17 -Wall -Wextra -g -DDEBUG
INCLUDES = -I.
LIBS = -lm

# Project name
PROJECT = genesis

# Source files
SOURCES = main.cpp genesis_core.cpp
HEADERS = genesis_core.h genesis_config.h

# Object files
OBJECTS = $(SOURCES:.cpp=.o)
DEBUG_OBJECTS = $(SOURCES:.cpp=_debug.o)

# Target executable
TARGET = $(PROJECT)
DEBUG_TARGET = $(PROJECT)_debug

# Default target
all: $(TARGET)

# Debug target
debug: $(DEBUG_TARGET)

# Release build
$(TARGET): $(OBJECTS)
	@echo "Linking release build..."
	$(CXX) $(CXXFLAGS) $(OBJECTS) -o $(TARGET) $(LIBS)
	@echo "Build complete: $(TARGET)"

# Debug build
$(DEBUG_TARGET): $(DEBUG_OBJECTS)
	@echo "Linking debug build..."
	$(CXX) $(DEBUG_FLAGS) $(DEBUG_OBJECTS) -o $(DEBUG_TARGET) $(LIBS)
	@echo "Debug build complete: $(DEBUG_TARGET)"

# Compile release object files
%.o: %.cpp $(HEADERS)
	@echo "Compiling $< (release)..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Compile debug object files
%_debug.o: %.cpp $(HEADERS)
	@echo "Compiling $< (debug)..."
	$(CXX) $(DEBUG_FLAGS) $(INCLUDES) -c $< -o $@

# Run the program
run: $(TARGET)
	@echo "Running $(TARGET)..."
	./$(TARGET)

# Run debug version
run-debug: $(DEBUG_TARGET)
	@echo "Running $(DEBUG_TARGET)..."
	./$(DEBUG_TARGET)

# Run quick test
quick: $(TARGET)
	@echo "Running quick test..."
	./$(TARGET) --quick

# Run batch parameter analysis
batch: $(TARGET)
	@echo "Running batch analysis..."
	./$(TARGET) --batch

# Memory check (requires valgrind)
memcheck: $(DEBUG_TARGET)
	@echo "Running memory check..."
	valgrind --leak-check=full --show-leak-kinds=all ./$(DEBUG_TARGET)

# Performance profiling (requires gprof)
profile: CXXFLAGS += -pg
profile: $(TARGET)
	@echo "Running performance profile..."
	./$(TARGET)
	gprof $(TARGET) gmon.out > profile_report.txt
	@echo "Profile report saved to profile_report.txt"

# Clean build files
clean:
	@echo "Cleaning build files..."
	rm -f $(OBJECTS) $(DEBUG_OBJECTS) $(TARGET) $(DEBUG_TARGET)
	rm -f *.o *_debug.o
	rm -f gmon.out profile_report.txt
	@echo "Clean complete"

# Install (copy to /usr/local/bin)
install: $(TARGET)
	@echo "Installing $(TARGET)..."
	cp $(TARGET) /usr/local/bin/
	@echo "Installation complete"

# Uninstall
uninstall:
	@echo "Uninstalling $(TARGET)..."
	rm -f /usr/local/bin/$(TARGET)
	@echo "Uninstall complete"

# Create project documentation
docs:
	@echo "Generating documentation..."
	@echo "# Genesis Project Documentation" > README.md
	@echo "" >> README.md
	@echo "## Building" >> README.md
	@echo "```bash" >> README.md
	@echo "make          # Release build" >> README.md
	@echo "make debug    # Debug build" >> README.md
	@echo "make run      # Build and run" >> README.md
	@echo "```" >> README.md
	@echo "Documentation generated: README.md"

# Package for distribution
package: clean
	@echo "Creating distribution package..."
	tar -czf $(PROJECT).tar.gz *.cpp *.h Makefile
	@echo "Package created: $(PROJECT).tar.gz"

# Check coding style (requires clang-format)
format:
	@echo "Formatting source code..."
	clang-format -i *.cpp *.h
	@echo "Code formatting complete"

# Static analysis (requires cppcheck)
analyze:
	@echo "Running static analysis..."
	cppcheck --enable=all --std=c++17 *.cpp *.h
	@echo "Static analysis complete"

# Help target
help:
	@echo "Genesis Project Makefile"
	@echo "========================"
	@echo ""
	@echo "Available targets:"
	@echo "  all       - Build release version (default)"
	@echo "  debug     - Build debug version"
	@echo "  run       - Build and run release version"
	@echo "  run-debug - Build and run debug version"
	@echo "  quick     - Run quick test (small parameters)"
	@echo "  batch     - Run batch parameter analysis"
	@echo "  memcheck  - Run memory leak check (requires valgrind)"
	@echo "  profile   - Generate performance profile (requires gprof)"
	@echo "  clean     - Remove all build files"
	@echo "  install   - Install to /usr/local/bin"
	@echo "  uninstall - Remove from /usr/local/bin"
	@echo "  docs      - Generate documentation"
	@echo "  package   - Create distribution package"
	@echo "  format    - Format source code (requires clang-format)"
	@echo "  analyze   - Run static analysis (requires cppcheck)"
	@echo "  help      - Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make                    # Build release version"
	@echo "  make debug run-debug    # Build and run debug version"
	@echo "  make clean all run      # Clean, build, and run"

# Phony targets
.PHONY: all debug run run-debug quick batch memcheck profile clean install uninstall docs package format analyze help

# Default goal
.DEFAULT_GOAL := all
