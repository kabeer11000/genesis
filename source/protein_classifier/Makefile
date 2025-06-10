# Protein Classifier Makefile

CXX = g++
CXXFLAGS = -std=c++17 -Wall -Wextra -O2 -Iinclude
SRCDIR = src
OBJDIR = obj
SOURCES = $(wildcard $(SRCDIR)/*.cpp)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
TARGET = protein_classifier

.PHONY: all clean run logs

all: $(TARGET)

$(TARGET): $(OBJECTS) | logs
	@echo "Linking $(TARGET)..."
	$(CXX) $(OBJECTS) -o $(TARGET)
	@echo "Build completed successfully!"

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

logs:
	mkdir -p logs

clean:
	@echo "Cleaning build files..."
	rm -rf $(OBJDIR) $(TARGET) logs/*.log

run: $(TARGET)
	@echo "Running protein classifier..."
	./$(TARGET)

debug: CXXFLAGS += -g -DDEBUG
debug: $(TARGET)

install: $(TARGET)
	cp $(TARGET) /usr/local/bin/

help:
	@echo "Available targets:"
	@echo "  all    - Build the project"
	@echo "  clean  - Remove build files"
	@echo "  run    - Build and run the program"
	@echo "  debug  - Build with debug symbols"
	@echo "  help   - Show this help message"
