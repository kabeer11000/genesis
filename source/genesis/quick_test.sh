#!/bin/bash
echo "Genesis Quick Test Script"
echo "========================="

echo "Building project..."
make clean
make

echo ""
echo "Running quick test..."
./genesis --quick

echo ""
echo "Quick test complete!"
