#ifndef CLASSIFIER_H
#define CLASSIFIER_H

#include "protein.h"
#include "database.h"
#include <vector>
#include <map>
#include <string>
#include <iomanip>

class ProteinClassifier {
private:
    std::vector<std::string> selectedFunctions;
    std::map<std::string, double> weights;
    double threshold;
    double accuracy;
    
public:
    ProteinClassifier();
    
    // Training
    void train(const ProteinDatabase& db);
    void evolveWeights(const ProteinDatabase& db, int generations = 100);
    
    // Classification
    std::string classify(const Protein& p);
    double computeScore(const Protein& p);
    
    // Evaluation
    double evaluateAccuracy(const ProteinDatabase& db);
    void testOnUnknownProteins(const ProteinDatabase& db);
    
    // Display
    void printClassifier() const;
    void printWeights() const;
};

#endif
