#ifndef PROTEIN_H
#define PROTEIN_H

#include <string>
#include <vector>
#include <iostream>
#include <map>

class Protein {
private:
    std::string sequence;
    std::string name;
    std::string classification;
    
public:
    Protein(const std::string& seq, const std::string& n, const std::string& cls = "unknown");
    
    // Basic getters
    std::string getSequence() const { return sequence; }
    std::string getName() const { return name; }
    std::string getClassification() const { return classification; }
    int getLength() const { return sequence.length(); }
    
    // Protein property calculations
    double calculateHydrophobicity() const;
    double calculateMolecularWeight() const;
    int countAminoAcid(char amino) const;
    double calculateChargeAtPH7() const;
    
    // Display functions
    void print() const;
    void printDetailed() const;
};

#endif
