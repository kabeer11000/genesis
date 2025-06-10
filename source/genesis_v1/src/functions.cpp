#include "../include/functions.h"
#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>

// Static member initialization
std::map<std::string, ProteinFunction> ProteinFunctions::fixedFunctions;
std::map<std::string, ProteinFunction> ProteinFunctions::dynamicFunctions;
int ProteinFunctions::dynamicFunctionCounter = 0;

void ProteinFunctions::initializeFixedFunctions() {
    std::cout << "[LOG] Initializing fixed functions..." << std::endl;
    
    fixedFunctions["f1"] = f1_hydrophobicity;
    fixedFunctions["f2"] = f2_length_normalized;
    fixedFunctions["f3"] = f3_charge_density;
    
    std::cout << "[LOG] Fixed functions initialized: f1, f2, f3" << std::endl;
}

double ProteinFunctions::f1_hydrophobicity(const Protein& p) {
    double hydro = p.calculateHydrophobicity();
    double result = std::abs(hydro) * 10.0; // Scale and normalize
    std::cout << "[LOG] f1_hydrophobicity(" << p.getName() << ") = " << result << std::endl;
    return result;
}

double ProteinFunctions::f2_length_normalized(const Protein& p) {
    double length = static_cast<double>(p.getLength());
    double result = (length / 100.0) * 20.0; // Normalize to 0-20 range
    std::cout << "[LOG] f2_length_normalized(" << p.getName() << ") = " << result << std::endl;
    return result;
}

double ProteinFunctions::f3_charge_density(const Protein& p) {
    double charge = std::abs(p.calculateChargeAtPH7());
    double length = static_cast<double>(p.getLength());
    double result = length > 0 ? (charge / length) * 100.0 : 0.0;
    std::cout << "[LOG] f3_charge_density(" << p.getName() << ") = " << result << std::endl;
    return result;
}

std::string ProteinFunctions::generateCombinationFunction(const std::string& func1, const std::string& func2, const std::string& operation) {
    std::string newFuncName = "dynamic_" + std::to_string(++dynamicFunctionCounter);
    
    std::cout << "[LOG] Generating combination function: " << newFuncName << " = " << func1 << " " << operation << " " << func2 << std::endl;
    
    ProteinFunction newFunc;
    
    if (operation == "+") {
        newFunc = [func1, func2](const Protein& p) -> double {
            return evaluateFunction(func1, p) + evaluateFunction(func2, p);
        };
    } else if (operation == "*") {
        newFunc = [func1, func2](const Protein& p) -> double {
            return evaluateFunction(func1, p) * evaluateFunction(func2, p) * 0.1; // Scale down
        };
    } else if (operation == "-") {
        newFunc = [func1, func2](const Protein& p) -> double {
            return std::abs(evaluateFunction(func1, p) - evaluateFunction(func2, p));
        };
    } else {
        newFunc = [func1, func2](const Protein& p) -> double {
            return (evaluateFunction(func1, p) + evaluateFunction(func2, p)) * 0.5;
        };
    }
    
    dynamicFunctions[newFuncName] = newFunc;
    std::cout << "[LOG] Dynamic function " << newFuncName << " created successfully" << std::endl;
    
    return newFuncName;
}

std::string ProteinFunctions::generateRandomFunction() {
    std::vector<std::string> baseFunc = {"f1", "f2", "f3"};
    std::vector<std::string> operations = {"+", "*", "-", "avg"};
    
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> funcDist(0, baseFunc.size() - 1);
    std::uniform_int_distribution<> opDist(0, operations.size() - 1);
    
    std::string func1 = baseFunc[funcDist(gen)];
    std::string func2 = baseFunc[funcDist(gen)];
    std::string op = operations[opDist(gen)];
    
    std::cout << "[LOG] Generating random function: " << func1 << " " << op << " " << func2 << std::endl;
    
    return generateCombinationFunction(func1, func2, op);
}

double ProteinFunctions::evaluateFunction(const std::string& funcName, const Protein& p) {
    // Check fixed functions first
    if (fixedFunctions.find(funcName) != fixedFunctions.end()) {
        return fixedFunctions[funcName](p);
    }
    
    // Check dynamic functions
    if (dynamicFunctions.find(funcName) != dynamicFunctions.end()) {
        return dynamicFunctions[funcName](p);
    }
    
    std::cout << "[LOG] Warning: Function " << funcName << " not found, returning 0.0" << std::endl;
    return 0.0;
}

std::vector<std::string> ProteinFunctions::getAllFunctionNames() {
    std::vector<std::string> names;
    
    for (const auto& pair : fixedFunctions) {
        names.push_back(pair.first);
    }
    
    for (const auto& pair : dynamicFunctions) {
        names.push_back(pair.first);
    }
    
    std::cout << "[LOG] Total functions available: " << names.size() << std::endl;
    return names;
}

std::map<std::string, double> ProteinFunctions::computeFunctionTable(const Protein& p) {
    std::map<std::string, double> table;
    
    std::cout << "[LOG] Computing function table for protein: " << p.getName() << std::endl;
    
    std::vector<std::string> allFunctions = getAllFunctionNames();
    
    for (const std::string& funcName : allFunctions) {
        double value = evaluateFunction(funcName, p);
        table[funcName] = value;
        std::cout << "[LOG] " << funcName << " = " << value << std::endl;
    }
    
    std::cout << "[LOG] Function table computed with " << table.size() << " entries" << std::endl;
    return table;
}

void ProteinFunctions::printFunctionTable(const Protein& p) {
    std::cout << "\n=== Function Table for " << p.getName() << " ===" << std::endl;
    std::map<std::string, double> table = computeFunctionTable(p);
    
    for (const auto& pair : table) {
        std::cout << std::setw(15) << pair.first << " : " << std::setw(10) << std::fixed << std::setprecision(3) << pair.second << std::endl;
    }
    std::cout << "==============================" << std::endl;
}
