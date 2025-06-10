#include "../include/database.h"
#include "../include/functions.h"
#include <iostream>
#include <iomanip>

void ProteinDatabase::initializeDatabase() {
    std::cout << "[LOG] Initializing protein database..." << std::endl;
    
    // Simple protein sequences (real amino acid sequences, simplified)
    proteins.clear();
    testProteins.clear();
    
    // Database proteins (4 simple small proteins)
    proteins.emplace_back("MKLLNVINFVFLMFVSSSRILGMENAMPAASLIQVVNTIIAKMKLDNN", "Insulin_A", "hormone");
    proteins.emplace_back("FVNQHLCGSHLVEALYLVCGERGFFYTPKA", "Insulin_B", "hormone");
    proteins.emplace_back("GSSGSSGKFGGGGGKFGGGKFNGRSPEPNNPKVTGGGKFGGGGKFGG", "Elastin", "structural");
    proteins.emplace_back("MAKKTSSAQKRGARRGWFSSRAKAKSSRSS", "Histone_H4", "regulatory");
    
    // Test protein (not in database)
    testProteins.emplace_back("MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKN", "Unknown_Test", "unknown");
    
    std::cout << "[LOG] Database initialized with " << proteins.size() << " known proteins" << std::endl;
    std::cout << "[LOG] Test set initialized with " << testProteins.size() << " unknown proteins" << std::endl;
}

void ProteinDatabase::addProtein(const Protein& p) {
    proteins.push_back(p);
    std::cout << "[LOG] Added protein " << p.getName() << " to database" << std::endl;
}

void ProteinDatabase::addTestProtein(const Protein& p) {
    testProteins.push_back(p);
    std::cout << "[LOG] Added protein " << p.getName() << " to test set" << std::endl;
}

void ProteinDatabase::printDatabase() const {
    std::cout << "\n=== Protein Database ===" << std::endl;
    std::cout << "Known Proteins (" << proteins.size() << "):" << std::endl;
    
    for (size_t i = 0; i < proteins.size(); ++i) {
        std::cout << "[" << i+1 << "] ";
        proteins[i].print();
        std::cout << std::endl;
    }
    
    std::cout << "\nTest Proteins (" << testProteins.size() << "):" << std::endl;
    for (size_t i = 0; i < testProteins.size(); ++i) {
        std::cout << "[T" << i+1 << "] ";
        testProteins[i].print();
        std::cout << std::endl;
    }
    std::cout << "========================" << std::endl;
}

void ProteinDatabase::computeAllTables() const {
    std::cout << "[LOG] Computing function tables for all proteins..." << std::endl;
    
    for (const Protein& p : proteins) {
        ProteinFunctions::printFunctionTable(p);
    }
    
    for (const Protein& p : testProteins) {
        ProteinFunctions::printFunctionTable(p);
    }
    
    std::cout << "[LOG] All function tables computed" << std::endl;
}

Protein* ProteinDatabase::findProtein(const std::string& name) {
    for (Protein& p : proteins) {
        if (p.getName() == name) {
            std::cout << "[LOG] Found protein: " << name << std::endl;
            return &p;
        }
    }
    
    for (Protein& p : testProteins) {
        if (p.getName() == name) {
            std::cout << "[LOG] Found test protein: " << name << std::endl;
            return &p;
        }
    }
    
    std::cout << "[LOG] Protein not found: " << name << std::endl;
    return nullptr;
}

void ProteinDatabase::printStatistics() const {
    std::cout << "\n=== Database Statistics ===" << std::endl;
    std::cout << "Total known proteins: " << proteins.size() << std::endl;
    std::cout << "Total test proteins: " << testProteins.size() << std::endl;
    
    if (!proteins.empty()) {
        double avgLength = 0.0;
        for (const Protein& p : proteins) {
            avgLength += p.getLength();
        }
        avgLength /= proteins.size();
        std::cout << "Average protein length: " << std::fixed << std::setprecision(1) << avgLength << std::endl;
    }
    
    std::cout << "===========================" << std::endl;
}
