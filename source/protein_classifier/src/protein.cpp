#include "../include/protein.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <map>

// Amino acid properties
static std::map<char, double> hydrophobicity = {
    {'A', 1.8}, {'R', -4.5}, {'N', -3.5}, {'D', -3.5}, {'C', 2.5},
    {'Q', -3.5}, {'E', -3.5}, {'G', -0.4}, {'H', -3.2}, {'I', 4.5},
    {'L', 3.8}, {'K', -3.9}, {'M', 1.9}, {'F', 2.8}, {'P', -1.6},
    {'S', -0.8}, {'T', -0.7}, {'W', -0.9}, {'Y', -1.3}, {'V', 4.2}
};

static std::map<char, double> molecularWeight = {
    {'A', 89.1}, {'R', 174.2}, {'N', 132.1}, {'D', 133.1}, {'C', 121.0},
    {'Q', 146.1}, {'E', 147.1}, {'G', 75.1}, {'H', 155.2}, {'I', 131.2},
    {'L', 131.2}, {'K', 146.2}, {'M', 149.2}, {'F', 165.2}, {'P', 115.1},
    {'S', 105.1}, {'T', 119.1}, {'W', 204.2}, {'Y', 181.2}, {'V', 117.1}
};

static std::map<char, double> charge = {
    {'A', 0}, {'R', 1}, {'N', 0}, {'D', -1}, {'C', 0},
    {'Q', 0}, {'E', -1}, {'G', 0}, {'H', 0.1}, {'I', 0},
    {'L', 0}, {'K', 1}, {'M', 0}, {'F', 0}, {'P', 0},
    {'S', 0}, {'T', 0}, {'W', 0}, {'Y', 0}, {'V', 0}
};

Protein::Protein(const std::string& seq, const std::string& n, const std::string& cls) 
    : sequence(seq), name(n), classification(cls) {
    std::cout << "[LOG] Created protein: " << name << " with sequence length: " << sequence.length() << std::endl;
}

double Protein::calculateHydrophobicity() const {
    double total = 0.0;
    int validCount = 0;
    
    std::cout << "[LOG] Calculating hydrophobicity for " << name << std::endl;
    
    for (char amino : sequence) {
        if (hydrophobicity.find(amino) != hydrophobicity.end()) {
            total += hydrophobicity[amino];
            validCount++;
        }
    }
    
    double result = validCount > 0 ? total / validCount : 0.0;
    std::cout << "[LOG] Hydrophobicity result: " << result << std::endl;
    return result;
}

double Protein::calculateMolecularWeight() const {
    double total = 0.0;
    
    std::cout << "[LOG] Calculating molecular weight for " << name << std::endl;
    
    for (char amino : sequence) {
        if (molecularWeight.find(amino) != molecularWeight.end()) {
            total += molecularWeight[amino];
        }
    }
    
    std::cout << "[LOG] Molecular weight result: " << total << " Da" << std::endl;
    return total;
}

int Protein::countAminoAcid(char amino) const {
    int count = std::count(sequence.begin(), sequence.end(), amino);
    std::cout << "[LOG] Count of " << amino << " in " << name << ": " << count << std::endl;
    return count;
}

double Protein::calculateChargeAtPH7() const {
    double total = 0.0;
    
    std::cout << "[LOG] Calculating charge at pH 7 for " << name << std::endl;
    
    for (char amino : sequence) {
        if (charge.find(amino) != charge.end()) {
            total += charge[amino];
        }
    }
    
    std::cout << "[LOG] Charge at pH 7 result: " << total << std::endl;
    return total;
}

void Protein::print() const {
    std::cout << "Protein: " << name << " [" << classification << "]" << std::endl;
    std::cout << "Sequence: " << sequence << std::endl;
    std::cout << "Length: " << getLength() << " amino acids" << std::endl;
}

void Protein::printDetailed() const {
    std::cout << "\n=== Detailed Protein Analysis ===" << std::endl;
    print();
    std::cout << "Hydrophobicity: " << calculateHydrophobicity() << std::endl;
    std::cout << "Molecular Weight: " << calculateMolecularWeight() << " Da" << std::endl;
    std::cout << "Charge at pH 7: " << calculateChargeAtPH7() << std::endl;
    std::cout << "=================================" << std::endl;
}
