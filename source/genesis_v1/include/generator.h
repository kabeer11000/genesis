#ifndef GENERATOR_H
#define GENERATOR_H

#include "functions.h"
#include "database.h"
#include <vector>
#include <string>

class FunctionGenerator {
private:
    static int populationSize;
    static std::vector<std::string> chromosome;
    
public:
    // Generator operations
    static void generatePopulation(int size = 200);
    static void evolveChromosome(const ProteinDatabase& db, int generations = 50);
    
    // Function creation
    static std::vector<std::string> generateRandomFunctions(int count);
    static std::string combineFunctions(const std::string& f1, const std::string& f2);
    
    // Evaluation
    static double evaluateFitness(const std::vector<std::string>& functions, const ProteinDatabase& db);
    static void selectBestFunctions(const ProteinDatabase& db);
    
    // Display
    static void printPopulation();
    static void printChromosome();
};

#endif
