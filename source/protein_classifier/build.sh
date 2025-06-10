#!/bin/bash

echo "=== Building Protein Classifier ==="

# Create necessary directories
mkdir -p obj logs data

# Build the project
echo "Compiling source files..."
make clean
make all

if [ $? -eq 0 ]; then
    echo "Build successful!"
    echo ""
    echo "To run the program:"
    echo "  make run"
    echo ""
    echo "To clean build files:"
    echo "  make clean"
    echo ""
    echo "Available commands:"
    echo "  ./protein_classifier    # Run directly"
    echo "  make debug             # Build with debug info"
    echo "  make help              # Show help"
else
    echo "Build failed!"
    exit 1
fi
