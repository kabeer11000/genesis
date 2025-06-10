#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "protein.h"
#include <functional>
#include <vector>
#include <map>
#include <string>

// Function type definition
typedef std::function<double(const Protein&)> ProteinFunction;

class ProteinFunctions {
private:
    static std::map<std::string, ProteinFunction> fixedFunctions;
    static std::map<std::string, ProteinFunction> dynamicFunctions;
    static int dynamicFunctionCounter;
    
public:
    // Initialize fixed functions
    static void initializeFixedFunctions();
    
    // Fixed function implementations
    static double f1_hydrophobicity(const Protein& p);
    static double f2_length_normalized(const Protein& p);
    static double f3_charge_density(const Protein& p);
    
    // Dynamic function generation
    static std::string generateCombinationFunction(const std::string& func1, const std::string& func2, const std::string& operation);
    static std::string generateRandomFunction();
    
    // Function evaluation
    static double evaluateFunction(const std::string& funcName, const Protein& p);
    static std::vector<std::string> getAllFunctionNames();
    
    // Function table operations
    static std::map<std::string, double> computeFunctionTable(const Protein& p);
    static void printFunctionTable(const Protein& p);
};

#endif
