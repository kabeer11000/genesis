#!/bin/bash

echo "=== Testing Protein Classifier ==="

# Build if necessary
if [ ! -f protein_classifier ]; then
    echo "Building project..."
    make all
fi

# Run the program and capture output
echo "Running protein classification..."
./protein_classifier > test_output.log 2>&1

# Check if execution was successful
if [ $? -eq 0 ]; then
    echo "Test execution successful!"
    echo ""
    echo "=== Key Results ==="
    echo "1. Database loaded:"
    grep "Database initialized" test_output.log
    
    echo ""
    echo "2. Functions generated:"
    grep "functions initialized" test_output.log
    grep "Dynamic function.*created" test_output.log | head -3
    
    echo ""
    echo "3. Classification results:"
    grep "Predicted class" test_output.log
    
    echo ""
    echo "=== Full Output ==="
    echo "Check test_output.log for complete execution log"
    echo "Check logs/execution.log for system logs"
else
    echo "Test execution failed!"
    echo "Check test_output.log for error details"
    exit 1
fi
