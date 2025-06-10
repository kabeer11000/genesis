#ifndef DATABASE_H
#define DATABASE_H

#include "protein.h"
#include <vector>
#include <string>
#include <map>

class ProteinDatabase {
private:
    std::vector<Protein> proteins;
    std::vector<Protein> testProteins;
    
public:
    // Database initialization
    void initializeDatabase();
    void addProtein(const Protein& p);
    void addTestProtein(const Protein& p);
    
    // Getters
    const std::vector<Protein>& getProteins() const { return proteins; }
    const std::vector<Protein>& getTestProteins() const { return testProteins; }
    
    // Database operations
    void printDatabase() const;
    void computeAllTables() const;
    Protein* findProtein(const std::string& name);
    
    // Statistics
    void printStatistics() const;
};

#endif
