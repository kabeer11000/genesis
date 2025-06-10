#ifndef GENESIS_CONFIG_H
#define GENESIS_CONFIG_H

#include "genesis_core.h"
#include <map>
#include <string>
#include <chrono>
#include <iomanip>

// Configuration structure for experiments
struct ExperimentConfig {
    // Evolutionary algorithm parameters
    size_t population_size = 100;
    double mutation_rate = 0.1;
    double mutation_strength = 0.5;
    double crossover_rate = 0.8;
    size_t generations = 100;
    
    // Function configuration
    bool use_angle_function = true;
    bool use_hydrophobicity_function = true;
    bool use_distance_constraint_function = true;
    bool use_combination_functions = true;
    
    // Test sequences
    std::vector<std::string> target_sequences;
    
    // Output settings
    bool verbose = true;
    bool save_results = false;
    std::string output_directory = "./results/";
    
    // Performance settings
    size_t max_threads = 1;
    bool use_gpu = false;
    
    ExperimentConfig() {
        // Default test sequences
        target_sequences = {
            "MKQLEDKVEELLSKNYHLENEVARLKKLVGER",  // 32 residues
            "AVILGLDKLKQKGDDDLKELDL",           // 21 residues
            "MVKSLQIGSCFGTVFGKY",              // 18 residues
            "GPLGS"                            // 5 residues (short test)
        };
    }
};

// Predefined experimental configurations
class ExperimentConfigs {
public:
    static ExperimentConfig getQuickTest() {
        ExperimentConfig config;
        config.population_size = 20;
        config.generations = 20;
        config.target_sequences = {"GPLGS", "MKQLED"};
        return config;
    }
    
    static ExperimentConfig getStandardTest() {
        ExperimentConfig config;
        config.population_size = 100;
        config.generations = 100;
        return config;
    }
    
    static ExperimentConfig getIntensiveTest() {
        ExperimentConfig config;
        config.population_size = 200;
        config.generations = 500;
        config.mutation_rate = 0.05;
        config.mutation_strength = 0.3;
        return config;
    }
    
    static ExperimentConfig getParameterSweep() {
        ExperimentConfig config;
        config.population_size = 50;
        config.generations = 50;
        config.target_sequences = {"MKQLEDKVEELLSKNYHLENEVARLKKLVGER"};
        return config;
    }
    
    static ExperimentConfig getMiniTest() {
        ExperimentConfig config;
        config.population_size = 10;
        config.generations = 10;
        config.target_sequences = {"GPLGS"};
        return config;
    }
};

// Experiment runner with configuration support
class ConfigurableExperiment {
private:
    ExperimentConfig config;
    std::unique_ptr<EvolutionaryAlgorithm> ea;
    std::unique_ptr<FunctionExecutor> executor;
    
public:
    ConfigurableExperiment(const ExperimentConfig& cfg) : config(cfg) {
        setupExperiment();
    }
    
    void setupExperiment();
    void run();
    void printConfig();
    void printResults(long long duration_ms);
    void saveResults();
    
    Chromosome getBestChromosome();
    std::vector<double> getPopulationFitness();
};

// Batch experiment runner for parameter sweeps
class BatchExperimentRunner {
public:
    static void runMutationRateSweep();
    static void runPopulationSizeSweep();
    static void runFunctionConfigurationTest();
};

#endif // GENESIS_CONFIG_H
