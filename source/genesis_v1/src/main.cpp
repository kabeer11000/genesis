#include <iomanip>
#include "../include/protein.h"
#include "../include/functions.h"
#include "../include/database.h"
#include "../include/classifier.h"
#include "../include/generator.h"
#include <iostream>
#include <fstream>


void logToFile(const std::string& message) {
    std::ofstream logFile("logs/execution.log", std::ios::app);
    logFile << message << std::endl;
    logFile.close();
}

int main() {
    std::cout << "=== Protein Classification System ===" << std::endl;
    std::cout << "Initializing..." << std::endl;
    
    // Create logs directory
    system("mkdir -p logs");
    
    // Initialize system
    ProteinFunctions::initializeFixedFunctions();
    
    // Create and initialize database
    ProteinDatabase database;
    database.initializeDatabase();
    
    std::cout << "\n1. Database Overview" << std::endl;
    database.printDatabase();
    database.printStatistics();
    
    std::cout << "\n2. Computing Function Tables" << std::endl;
    database.computeAllTables();
    
    std::cout << "\n3. Generating Function Population" << std::endl;
    FunctionGenerator::generatePopulation(50); // Smaller for demo
    FunctionGenerator::printPopulation();
    
    std::cout << "\n4. Evolving Functions" << std::endl;
    FunctionGenerator::evolveChromosome(database, 20);
    FunctionGenerator::selectBestFunctions(database);
    
    std::cout << "\n5. Training Classifier" << std::endl;
    ProteinClassifier classifier;
    classifier.train(database);
    classifier.printClassifier();
    classifier.printWeights();
    
    std::cout << "\n6. Testing on Unknown Proteins" << std::endl;
    classifier.testOnUnknownProteins(database);
    
    std::cout << "\n7. Final Analysis" << std::endl;
    std::cout << "System successfully classified proteins using:" << std::endl;
    std::cout << "- Fixed functions (f1, f2, f3)" << std::endl;
    std::cout << "- Dynamically generated functions" << std::endl;
    std::cout << "- Evolutionary optimization" << std::endl;
    std::cout << "- Function table intermediate representation" << std::endl;
    
    // Log completion
    logToFile("Protein classification system execution completed successfully");
    
    std::cout << "\n=== Execution Complete ===" << std::endl;
    std::cout << "Check logs/execution.log for detailed execution log" << std::endl;
    
    return 0;
}
