#include "genesis_core.h"
#include "genesis_config.h"
#include <iostream>
#include <chrono>
#include <iomanip>

void runBasicEvolution() {
    std::cout << "=== Genesis Protein Evolution Experiment ===" << std::endl;
    
    // Create evolutionary algorithm
    size_t population_size = 100;
    double mutation_rate = 0.1;
    double mutation_strength = 0.5;
    double crossover_rate = 0.8;
    
    EvolutionaryAlgorithm ea(population_size, mutation_rate, mutation_strength, crossover_rate);
    
    // Initialize with the total weight count needed for all functions
    FunctionExecutor executor = FunctionExecutor::createStandardSet();
    size_t weight_count = executor.getTotalWeightCount();
    std::cout << "Total weights in chromosome: " << weight_count << std::endl;
    
    ea.initialize(weight_count);
    
    // Define target protein sequences to evolve towards
    std::vector<std::string> target_sequences = {
        "MKQLEDKVEELLSKNYHLENEVARLKKLVGER",  // Small protein
        "AVILGLDKLKQKGDDDLKELDL",           // Another test sequence
        "MVKSLQIGSCFGTVFGKY"                // Short peptide
    };
    
    std::cout << "Target sequences:" << std::endl;
    for (size_t i = 0; i < target_sequences.size(); ++i) {
        std::cout << "  " << i+1 << ": " << target_sequences[i] << " (length: " 
                  << target_sequences[i].length() << ")" << std::endl;
    }
    
    // Time the evolution
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Run evolution
    size_t generations = 100;
    std::cout << "\nStarting evolution for " << generations << " generations..." << std::endl;
    ea.evolve(generations, target_sequences);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    // Get best chromosome
    Chromosome best = ea.getBest();
    std::cout << "\n=== Results ===" << std::endl;
    std::cout << "Best fitness: " << std::fixed << std::setprecision(6) << best.fitness << std::endl;
    std::cout << "Evolution time: " << duration.count() << " ms" << std::endl;
    
    // Show population fitness statistics
    std::vector<double> fitness_values = ea.getPopulationFitness();
    double sum_fitness = 0.0;
    double max_fitness = 0.0;
    double min_fitness = 1000.0;
    
    for (double f : fitness_values) {
        sum_fitness += f;
        max_fitness = std::max(max_fitness, f);
        min_fitness = std::min(min_fitness, f);
    }
    double avg_fitness = sum_fitness / fitness_values.size();
    
    std::cout << "Population statistics:" << std::endl;
    std::cout << "  Average fitness: " << std::fixed << std::setprecision(6) << avg_fitness << std::endl;
    std::cout << "  Max fitness: " << max_fitness << std::endl;
    std::cout << "  Min fitness: " << min_fitness << std::endl;
    
    // Test the best chromosome on each target sequence
    std::cout << "\n=== Testing Best Chromosome ===" << std::endl;
    for (size_t i = 0; i < target_sequences.size(); ++i) {
        Protein input(target_sequences[i]);
        Protein predicted = executor.execute(input, best);
        
        std::cout << "Sequence " << i+1 << ": " << target_sequences[i] << std::endl;
        std::cout << "  Predicted structure coordinates (first 5 residues):" << std::endl;
        
        for (size_t j = 0; j < std::min(5UL, predicted.backbone_coords.size()); ++j) {
            const auto& coord = predicted.backbone_coords[j];
            std::cout << "    " << j+1 << ": (" 
                      << std::fixed << std::setprecision(2) 
                      << coord.x << ", " << coord.y << ", " << coord.z << ")" << std::endl;
        }
        
        // Show some features
        if (predicted.phi_angles.size() > 0) {
            std::cout << "  Phi angle (residue 2): " 
                      << std::fixed << std::setprecision(2) 
                      << predicted.phi_angles[1] * 180.0 / M_PI << "°" << std::endl;
        }
        std::cout << std::endl;
    }
    
    // Show best chromosome weights (first 10)
    std::cout << "Best chromosome weights (first 10):" << std::endl;
    for (size_t i = 0; i < std::min(10UL, best.weights.size()); ++i) {
        std::cout << "  w[" << i << "] = " << std::fixed << std::setprecision(4) 
                  << best.weights[i] << std::endl;
    }
}

void runParameterTesting() {
    std::cout << "\n=== Parameter Testing ===" << std::endl;
    
    std::vector<double> mutation_rates = {0.05, 0.1, 0.2};
    std::vector<size_t> population_sizes = {50, 100, 200};
    
    std::string test_sequence = "MKQLEDKVEELLSKNYHLENEVARLKKLVGER";
    std::vector<std::string> target_sequences = {test_sequence};
    
    std::cout << "Testing different parameters on sequence: " << test_sequence << std::endl;
    std::cout << std::setw(10) << "PopSize" << std::setw(10) << "MutRate" 
              << std::setw(15) << "BestFitness" << std::setw(10) << "Time(ms)" << std::endl;
    std::cout << std::string(45, '-') << std::endl;
    
    for (size_t pop_size : population_sizes) {
        for (double mut_rate : mutation_rates) {
            EvolutionaryAlgorithm ea(pop_size, mut_rate, 0.5, 0.8);
            
            FunctionExecutor executor = FunctionExecutor::createStandardSet();
            size_t weight_count = executor.getTotalWeightCount();
            ea.initialize(weight_count);
            
            auto start_time = std::chrono::high_resolution_clock::now();
            ea.evolve(50, target_sequences); // Shorter evolution for testing
            auto end_time = std::chrono::high_resolution_clock::now();
            
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
            Chromosome best = ea.getBest();
            
            std::cout << std::setw(10) << pop_size 
                      << std::setw(10) << std::fixed << std::setprecision(2) << mut_rate
                      << std::setw(15) << std::setprecision(6) << best.fitness
                      << std::setw(10) << duration.count() << std::endl;
        }
    }
}

void demonstrateGraphRepresentation() {
    std::cout << "\n=== Function Graph Representation ===" << std::endl;
    
    FunctionExecutor executor = FunctionExecutor::createStandardSet();
    
    std::cout << "Function dependency graph:" << std::endl;
    std::cout << "Node 0: AngleFunction (inputs: protein_state) -> phi/psi angles" << std::endl;
    std::cout << "Node 1: HydrophobicityFunction (inputs: protein_state) -> hydrophobic positioning" << std::endl;
    std::cout << "Node 2: DistanceConstraintFunction (inputs: protein_state) -> distance corrections" << std::endl;
    std::cout << "Node 3: CombinationFunction(0,1) (inputs: Node0, Node1) -> combined angle/hydro" << std::endl;
    std::cout << "Node 4: CombinationFunction(1,2) (inputs: Node1, Node2) -> combined hydro/distance" << std::endl;
    std::cout << "Node 5: CombinationFunction(0,2,3) (inputs: Node0, Node2, Node3) -> final combination" << std::endl;
    for (int i = 6; i < N; i++) {
        executor.
    }
    
    std::cout << "\nData flow: Input -> [0,1,2] -> [3,4] -> [5] -> Output" << std::endl;
    std::cout << "Total weights needed: " << executor.getTotalWeightCount() << std::endl;
}

void runConfigurableExperiments() {
    std::cout << "\n=== Configurable Experiments ===" << std::endl;
    
    // Quick test
    std::cout << "Running quick test..." << std::endl;
    ConfigurableExperiment quick_exp(ExperimentConfigs::getQuickTest());
    quick_exp.run();
    
    // Standard test
    std::cout << "Running standard test..." << std::endl;
    ConfigurableExperiment standard_exp(ExperimentConfigs::getStandardTest());
    standard_exp.run();
}

void runBatchExperiments() {
    std::cout << "\n=== Batch Parameter Analysis ===" << std::endl;
    
    // Run parameter sweeps
    BatchExperimentRunner::runMutationRateSweep();
    BatchExperimentRunner::runPopulationSizeSweep();
    BatchExperimentRunner::runFunctionConfigurationTest();
}

int main(int argc, char* argv[]) {
    std::cout << "Genesis Protein Structure Prediction using Evolutionary Algorithms" << std::endl;
    std::cout << "=================================================================" << std::endl;
    
    try {
        // Check command line arguments for different modes
        bool run_quick = false;
        bool run_batch = false;
        bool run_all = true;
        
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "--quick" || arg == "-q") {
                run_quick = true;
                run_all = false;
            } else if (arg == "--batch" || arg == "-b") {
                run_batch = true;
                run_all = false;
            } else if (arg == "--help" || arg == "-h") {
                std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
                std::cout << "Options:" << std::endl;
                std::cout << "  --quick, -q    Run quick test only" << std::endl;
                std::cout << "  --batch, -b    Run batch parameter analysis" << std::endl;
                std::cout << "  --help, -h     Show this help message" << std::endl;
                std::cout << "  (no options)   Run full demonstration" << std::endl;
                return 0;
            }
        }
        
        if (run_quick) {
            std::cout << "Running quick test mode..." << std::endl;
            ConfigurableExperiment quick_exp(ExperimentConfigs::getQuickTest());
            quick_exp.run();
        } else if (run_batch) {
            runBatchExperiments();
        } else if (run_all) {
            // Main evolution experiment
            runBasicEvolution();
            
            // Parameter testing
            runParameterTesting();
            
            // Configurable experiments
            runConfigurableExperiments();
            
            // Batch experiments (shortened for demo)
            std::cout << "\n=== Sample Batch Analysis ===" << std::endl;
            BatchExperimentRunner::runFunctionConfigurationTest();
            
            // Show graph representation
            demonstrateGraphRepresentation();
        }
        
        std::cout << "\n=== Experiment Complete ===" << std::endl;
        std::cout << "The evolutionary algorithm has successfully evolved function weights" << std::endl;
        std::cout << "to predict protein structures. You can now:" << std::endl;
        std::cout << "1. Add more sophisticated functions" << std::endl;
        std::cout << "2. Load real PDB files for training" << std::endl;
        std::cout << "3. Experiment with different evolutionary parameters" << std::endl;
        std::cout << "4. Add neural networks or PSO algorithms" << std::endl;
        std::cout << "5. Use configuration system for systematic experiments" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}
