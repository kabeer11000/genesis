#include "../include/generator.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <iomanip>

// Static member initialization
int FunctionGenerator::populationSize = 200;
std::vector<std::string> FunctionGenerator::chromosome;

void FunctionGenerator::generatePopulation(int size) {
    populationSize = size;
    chromosome.clear();
    
    std::cout << "[LOG] Generating function population of size " << size << std::endl;
    
    // Start with base functions
    chromosome.push_back("f1");
    chromosome.push_back("f2");
    chromosome.push_back("f3");
    
    // Generate random combinations
    for (int i = 3; i < size; ++i) {
        std::string newFunc = ProteinFunctions::generateRandomFunction();
        chromosome.push_back(newFunc);
        
        if (i % 50 == 0) {
            std::cout << "[LOG] Generated " << i << " functions..." << std::endl;
        }
    }
    
    std::cout << "[LOG] Population generation completed with " << chromosome.size() << " functions" << std::endl;
}

void FunctionGenerator::evolveChromosome(const ProteinDatabase& db, int generations) {
    std::cout << "[LOG] Evolving chromosome over " << generations << " generations..." << std::endl;
    
    std::random_device rd;
    std::mt19937 rng(rd());
    
    for (int generation = 0; generation < generations; ++generation) {
        // Create offspring by combining existing functions
        std::vector<std::string> offspring;
        
        std::uniform_int_distribution<> dis(0, chromosome.size() - 1);
        
        for (int i = 0; i < 10; ++i) { // Generate 10 offspring per generation
            int idx1 = dis(rng);
            int idx2 = dis(rng);
            
            if (idx1 != idx2) {
                std::string newFunc = combineFunctions(chromosome[idx1], chromosome[idx2]);
                offspring.push_back(newFunc);
            }
        }
        
        // Add offspring to chromosome
        chromosome.insert(chromosome.end(), offspring.begin(), offspring.end());
        
        // Selection: keep best functions (simplified - keep recent ones)
        if (chromosome.size() > populationSize * 1.2) {
            chromosome.resize(populationSize);
        }
        
        if (generation % 10 == 0) {
            std::cout << "[LOG] Generation " << generation << ": Population size = " << chromosome.size() << std::endl;
        }
    }
    
    std::cout << "[LOG] Chromosome evolution completed" << std::endl;
}

std::vector<std::string> FunctionGenerator::generateRandomFunctions(int count) {
    std::vector<std::string> functions;
    
    std::cout << "[LOG] Generating " << count << " random functions..." << std::endl;
    
    for (int i = 0; i < count; ++i) {
        std::string func = ProteinFunctions::generateRandomFunction();
        functions.push_back(func);
    }
    
    std::cout << "[LOG] Random function generation completed" << std::endl;
    return functions;
}

std::string FunctionGenerator::combineFunctions(const std::string& f1, const std::string& f2) {
    std::vector<std::string> operations = {"+", "*", "-", "avg"};
    
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_int_distribution<> dis(0, operations.size() - 1);
    
    std::string op = operations[dis(rng)];
    
    std::cout << "[LOG] Combining functions: " << f1 << " " << op << " " << f2 << std::endl;
    
    return ProteinFunctions::generateCombinationFunction(f1, f2, op);
}

double FunctionGenerator::evaluateFitness(const std::vector<std::string>& functions, const ProteinDatabase& db) {
    // Simple fitness: variance of function values across proteins
    double fitness = 0.0;
    
    const std::vector<Protein>& proteins = db.getProteins();
    
    for (const std::string& func : functions) {
        std::vector<double> values;
        
        for (const Protein& p : proteins) {
            double val = ProteinFunctions::evaluateFunction(func, p);
            values.push_back(val);
        }
        
        // Calculate variance
        if (!values.empty()) {
            double mean = 0.0;
            for (double v : values) mean += v;
            mean /= values.size();
            
            double variance = 0.0;
            for (double v : values) {
                variance += (v - mean) * (v - mean);
            }
            variance /= values.size();
            
            fitness += variance;
        }
    }
    
    std::cout << "[LOG] Fitness evaluation result: " << fitness << std::endl;
    return fitness;
}

void FunctionGenerator::selectBestFunctions(const ProteinDatabase& db) {
    std::cout << "[LOG] Selecting best functions from population..." << std::endl;
    
    // Simple selection: keep functions with good discriminative power
    std::vector<std::pair<double, std::string>> functionScores;
    
    for (const std::string& func : chromosome) {
        std::vector<std::string> singleFunc = {func};
        double score = evaluateFitness(singleFunc, db);
        functionScores.push_back({score, func});
    }
    
    // Sort by score (higher is better)
    std::sort(functionScores.begin(), functionScores.end(), std::greater<>());
    
    // Keep top 50%
    chromosome.clear();
    for (size_t i = 0; i < functionScores.size() / 2; ++i) {
        chromosome.push_back(functionScores[i].second);
    }
    
    std::cout << "[LOG] Selected " << chromosome.size() << " best functions" << std::endl;
}

void FunctionGenerator::printPopulation() {
    std::cout << "\n=== Function Population ===" << std::endl;
    std::cout << "Population size: " << chromosome.size() << std::endl;
    
    for (size_t i = 0; i < std::min(chromosome.size(), static_cast<size_t>(10)); ++i) {
        std::cout << "[" << i+1 << "] " << chromosome[i] << std::endl;
    }
    
    if (chromosome.size() > 10) {
        std::cout << "... and " << (chromosome.size() - 10) << " more functions" << std::endl;
    }
    
    std::cout << "===========================" << std::endl;
}

void FunctionGenerator::printChromosome() {
    std::cout << "\n=== Current Chromosome ===" << std::endl;
    for (size_t i = 0; i < chromosome.size(); ++i) {
        std::cout << chromosome[i];
        if (i < chromosome.size() - 1) std::cout << ", ";
        if ((i + 1) % 5 == 0) std::cout << std::endl;
    }
    std::cout << "\n===========================" << std::endl;
}
