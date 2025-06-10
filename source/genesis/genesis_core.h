#ifndef GENESIS_CORE_H
#define GENESIS_CORE_H

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <random>
#include <fstream>
#include <cmath>

// Forward declarations
class Protein;
class Function;
class Chromosome;

// 3D coordinate structure
struct Coord3D {
    double x, y, z;
    Coord3D(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}
    
    Coord3D operator+(const Coord3D& other) const {
        return Coord3D(x + other.x, y + other.y, z + other.z);
    }
    
    Coord3D operator*(double scalar) const {
        return Coord3D(x * scalar, y * scalar, z * scalar);
    }
    
    double distance(const Coord3D& other) const {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return sqrt(dx*dx + dy*dy + dz*dz);
    }
};

// Amino acid properties
struct AminoAcid {
    char code;           // Single letter code
    double hydrophobicity;
    double volume;
    double charge;
    double mass;
    
    AminoAcid(char c = 'A');
};

// Protein representation
class Protein {
public:
    std::string sequence;
    std::vector<Coord3D> backbone_coords;  // CA coordinates
    std::vector<Coord3D> all_atom_coords;  // All heavy atoms
    std::vector<AminoAcid> amino_acids;
    
    // Structural features (computed by functions)
    std::vector<double> phi_angles;
    std::vector<double> psi_angles;
    std::vector<double> local_features;    // For storing intermediate results
    
    Protein(const std::string& seq);
    
    // Load from PDB file
    static Protein fromPDB(const std::string& filename);
    
    // Calculate RMSD between two protein structures
    double calculateRMSD(const Protein& other) const;
    
    // Get amino acid at position
    const AminoAcid& getAminoAcid(size_t pos) const;
    
    // Update coordinates
    void updateCoordinates(const std::vector<Coord3D>& coords);
    
    size_t length() const { return sequence.length(); }
};

// Base class for all functions
class Function {
public:
    virtual ~Function() = default;
    virtual void apply(Protein& protein, const std::vector<double>& weights) = 0;
    virtual std::string getName() const = 0;
    virtual size_t getWeightCount() const = 0;
};

// Predefined function: Calculates phi/psi angles
class AngleFunction : public Function {
public:
    void apply(Protein& protein, const std::vector<double>& weights) override;
    std::string getName() const override { return "AngleFunction"; }
    size_t getWeightCount() const override { return 2; } // weight for phi and psi
};

// Predefined function: Hydrophobicity-based positioning
class HydrophobicityFunction : public Function {
public:
    void apply(Protein& protein, const std::vector<double>& weights) override;
    std::string getName() const override { return "HydrophobicityFunction"; }
    size_t getWeightCount() const override { return 3; } // x, y, z influence
};

// Dynamic function: Linear combination of previous functions
class CombinationFunction : public Function {
private:
    std::vector<size_t> input_function_indices;
    
public:
    CombinationFunction(const std::vector<size_t>& indices) 
        : input_function_indices(indices) {}
    
    void apply(Protein& protein, const std::vector<double>& weights) override;
    std::string getName() const override { return "CombinationFunction"; }
    size_t getWeightCount() const override { return input_function_indices.size(); }
};

// Distance constraint function
class DistanceConstraintFunction : public Function {
public:
    void apply(Protein& protein, const std::vector<double>& weights) override;
    std::string getName() const override { return "DistanceConstraintFunction"; }
    size_t getWeightCount() const override { return 4; } // min_dist, max_dist, strength, range
};

// Chromosome represents weights for all functions
class Chromosome {
public:
    std::vector<double> weights;
    double fitness;
    
    Chromosome(size_t weight_count);
    Chromosome(const std::vector<double>& w) : weights(w), fitness(0.0) {}
    
    // Genetic operations
    Chromosome crossover(const Chromosome& other) const;
    void mutate(double mutation_rate, double mutation_strength);
    
    // Comparison for sorting
    bool operator<(const Chromosome& other) const {
        return fitness > other.fitness; // Higher fitness is better
    }
};

// Function executor - applies all functions sequentially
class FunctionExecutor {
private:
    std::vector<std::unique_ptr<Function>> functions;
    std::vector<size_t> weight_offsets; // Starting index for each function's weights
    
public:
    void addFunction(std::unique_ptr<Function> func);
    Protein execute(const Protein& input, const Chromosome& chromosome);
    size_t getTotalWeightCount() const;
    
    // Create a standard function set
    static FunctionExecutor createStandardSet();
};

// Fitness evaluator
class FitnessEvaluator {
private:
    std::vector<Protein> reference_proteins;
    
public:
    void addReference(const Protein& protein);
    void loadReferencesFromPDBDirectory(const std::string& directory);
    
    double evaluate(const Protein& predicted, const std::string& target_sequence);
    double evaluateStructuralFeatures(const Protein& predicted);
};

// Main evolutionary algorithm
class EvolutionaryAlgorithm {
private:
    size_t population_size;
    double mutation_rate;
    double mutation_strength;
    double crossover_rate;
    
    std::vector<Chromosome> population;
    FunctionExecutor executor;
    FitnessEvaluator evaluator;
    std::mt19937 rng;
    
public:
    EvolutionaryAlgorithm(size_t pop_size, double mut_rate, double mut_strength, double cross_rate);
    
    void initialize(size_t weight_count);
    void evolve(size_t generations, const std::vector<std::string>& target_sequences);
    
    Chromosome getBest() const;
    std::vector<double> getPopulationFitness() const;
    
    // Configuration
    void setMutationRate(double rate) { mutation_rate = rate; }
    void setMutationStrength(double strength) { mutation_strength = strength; }
    void setCrossoverRate(double rate) { crossover_rate = rate; }
    
    // Add reference proteins for fitness evaluation
    void addReferenceProtein(const Protein& protein);
    void loadReferencePDBs(const std::string& directory);
};

#endif // GENESIS_CORE_H
