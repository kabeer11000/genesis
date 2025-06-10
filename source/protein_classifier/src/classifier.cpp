#include "../include/classifier.h"
#include "../include/functions.h"
#include <iostream>
#include <random>
#include <algorithm>
#include <cmath>
#include <iomanip>

ProteinClassifier::ProteinClassifier() : threshold(0.5), accuracy(0.0) {
    std::cout << "[LOG] Protein classifier initialized" << std::endl;
}

void ProteinClassifier::train(const ProteinDatabase& db) {
    std::cout << "[LOG] Training classifier on database..." << std::endl;
    
    // Select functions to use
    selectedFunctions = {"f1", "f2", "f3"};
    
    // Generate some dynamic functions
    selectedFunctions.push_back(ProteinFunctions::generateRandomFunction());
    selectedFunctions.push_back(ProteinFunctions::generateRandomFunction());
    
    std::cout << "[LOG] Selected " << selectedFunctions.size() << " functions for classification" << std::endl;
    
    // Initialize weights randomly
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-1.0, 1.0);
    
    for (const std::string& func : selectedFunctions) {
        weights[func] = dis(gen);
        std::cout << "[LOG] Initial weight for " << func << ": " << weights[func] << std::endl;
    }
    
    // Evolve weights
    evolveWeights(db, 50);
    
    std::cout << "[LOG] Classifier training completed" << std::endl;
}

void ProteinClassifier::evolveWeights(const ProteinDatabase& db, int generations) {
    std::cout << "[LOG] Evolving weights over " << generations << " generations..." << std::endl;
    
    std::random_device rd;
    std::mt19937 rng(rd());
    std::uniform_real_distribution<> dis(-0.1, 0.1);
    
    double bestAccuracy = evaluateAccuracy(db);
    std::map<std::string, double> bestWeights = weights;
    
    for (int generation = 0; generation < generations; ++generation) {
        // Mutate weights
        std::map<std::string, double> newWeights = weights;
        for (auto& pair : newWeights) {
            pair.second += dis(rng);
        }
        
        // Test new weights
        std::map<std::string, double> oldWeights = weights;
        weights = newWeights;
        double newAccuracy = evaluateAccuracy(db);
        
        if (newAccuracy > bestAccuracy) {
            bestAccuracy = newAccuracy;
            bestWeights = newWeights;
            std::cout << "[LOG] Generation " << generation << ": New best accuracy = " << bestAccuracy << std::endl;
        } else {
            weights = oldWeights; // Revert
        }
        
        if (generation % 10 == 0) {
            std::cout << "[LOG] Generation " << generation << ": Current accuracy = " << evaluateAccuracy(db) << std::endl;
        }
    }
    
    weights = bestWeights;
    accuracy = bestAccuracy;
    
    std::cout << "[LOG] Evolution completed. Final accuracy: " << accuracy << std::endl;
}

std::string ProteinClassifier::classify(const Protein& p) {
    double score = computeScore(p);
    std::cout << "[LOG] Classification score for " << p.getName() << ": " << score << std::endl;
    
    if (score > threshold) {
        return "structural";
    } else {
        return "hormone";
    }
}

double ProteinClassifier::computeScore(const Protein& p) {
    double score = 0.0;
    
    for (const std::string& func : selectedFunctions) {
        double funcValue = ProteinFunctions::evaluateFunction(func, p);
        score += weights[func] * funcValue;
    }
    
    // Normalize score to 0-1 range
    score = 1.0 / (1.0 + std::exp(-score / 10.0));
    
    return score;
}

double ProteinClassifier::evaluateAccuracy(const ProteinDatabase& db) {
    const std::vector<Protein>& proteins = db.getProteins();
    int correct = 0;
    int total = proteins.size();
    
    for (const Protein& p : proteins) {
        std::string predicted = classify(p);
        std::string actual = p.getClassification();
        
        if (predicted == actual) {
            correct++;
        }
    }
    
    double acc = total > 0 ? static_cast<double>(correct) / total : 0.0;
    std::cout << "[LOG] Accuracy evaluation: " << correct << "/" << total << " = " << acc << std::endl;
    
    return acc;
}

void ProteinClassifier::testOnUnknownProteins(const ProteinDatabase& db) {
    std::cout << "\n=== Testing on Unknown Proteins ===" << std::endl;
    
    const std::vector<Protein>& testProteins = db.getTestProteins();
    
    for (const Protein& p : testProteins) {
        std::string prediction = classify(p);
        double score = computeScore(p);
        
        std::cout << "Protein: " << p.getName() << std::endl;
        std::cout << "Predicted class: " << prediction << std::endl;
        std::cout << "Confidence score: " << std::fixed << std::setprecision(3) << score << std::endl;
        std::cout << "---" << std::endl;
    }
    
    std::cout << "===================================" << std::endl;
}

void ProteinClassifier::printClassifier() const {
    std::cout << "\n=== Classifier Configuration ===" << std::endl;
    std::cout << "Threshold: " << threshold << std::endl;
    std::cout << "Accuracy: " << accuracy << std::endl;
    std::cout << "Selected functions: ";
    for (const std::string& func : selectedFunctions) {
        std::cout << func << " ";
    }
    std::cout << std::endl;
    std::cout << "================================" << std::endl;
}

void ProteinClassifier::printWeights() const {
    std::cout << "\n=== Classifier Weights ===" << std::endl;
    for (const auto& pair : weights) {
        std::cout << std::setw(15) << pair.first << " : " << std::setw(8) << std::fixed << std::setprecision(4) << pair.second << std::endl;
    }
    std::cout << "===========================" << std::endl;
}
