# Genesis: Evolutionary Protein Structure Prediction

A modular framework for evolving protein structure prediction functions using evolutionary algorithms, with planned support for neural networks and particle swarm optimization.

## Project Overview

Genesis implements a novel approach to protein structure prediction by evolving a composition of mathematical functions that sequentially transform amino acid sequences into 3D structural predictions. The system uses:

- **Sequential Function Application**: Functions are applied in sequence: `F(protein) = f_n(f_{n-1}(...f_1(f_0(protein))...))`
- **Evolutionary Optimization**: Genetic algorithms evolve the weights for each function to maximize prediction accuracy
- **Modular Architecture**: Easy to add new functions, algorithms, and evaluation metrics
- **Graph-based Representation**: Functions can depend on outputs from previous functions

## Quick Start

```bash
# Clone or extract the project
cd genesis/

# Build the project
make

# Run full demonstration
make run

# Quick test (small parameters)
make quick

# Parameter analysis
make batch

# See all options
make help
```

## Architecture

### Core Components

1. **Protein Class**: Represents protein sequences and 3D structures
2. **Function Classes**: Base class for all transformation functions
3. **Chromosome Class**: Represents evolved parameter vectors
4. **FunctionExecutor**: Manages sequential function application
5. **EvolutionaryAlgorithm**: Implements genetic algorithm optimization
6. **FitnessEvaluator**: Compares predicted vs. known structures

### Function Types

#### Predefined Functions
- **AngleFunction**: Calculates and applies phi/psi dihedral angles
- **HydrophobicityFunction**: Positions residues based on hydrophobicity
- **DistanceConstraintFunction**: Enforces realistic inter-residue distances

#### Dynamic Functions
- **CombinationFunction**: Linear combinations of previous function outputs
- **Custom Functions**: Easy to add domain-specific transformations

### Mathematical Representation

Each function can be represented as:
```
f_i(protein_state, weights) = modified_protein_state
```

The complete system represents:
```
F(sequence) = f_n(w_n, f_{n-1}(w_{n-1}, ... f_1(w_1, f_0(w_0, sequence))))
```

Where `w_i` are the evolved weights for function `f_i`.

## Usage Examples

### Basic Evolution
```cpp
#include "genesis_core.h"

// Create evolutionary algorithm
EvolutionaryAlgorithm ea(population_size=100, mutation_rate=0.1, 
                        mutation_strength=0.5, crossover_rate=0.8);

// Define target sequences
std::vector<std::string> targets = {"MKQLEDKVEELLSKNYHLENEVARLKKLVGER"};

// Initialize and evolve
FunctionExecutor executor = FunctionExecutor::createStandardSet();
ea.initialize(executor.getTotalWeightCount());
ea.evolve(generations=100, targets);

// Get best solution
Chromosome best = ea.getBest();
```

### Adding Custom Functions
```cpp
class CustomFunction : public Function {
public:
    void apply(Protein& protein, const std::vector<double>& weights) override {
        // Your custom transformation logic here
    }
    
    std::string getName() const override { return "CustomFunction"; }
    size_t getWeightCount() const override { return 3; }
};

// Add to executor
executor.addFunction(std::make_unique<CustomFunction>());
```

## Command Line Options

```bash
./genesis                 # Full demonstration
./genesis --quick         # Quick test (small parameters)
./genesis --batch         # Parameter analysis
./genesis --help          # Show help
```

## Project Structure

```
genesis/
├── genesis_core.h          # Header with class definitions
├── genesis_core.cpp        # Implementation file
├── genesis_config.h        # Configuration system
├── main.cpp               # Example usage and testing
├── Makefile              # Build configuration
├── README.md             # This file
├── data/                 # PDB files for reference (optional)
│   └── pdb_files/
└── results/              # Output directory
```

## Building Requirements

- C++17 compatible compiler (GCC 7+ or Clang 5+)
- Make build system
- Standard C++ libraries

## Next Steps

1. **Add Neural Networks**: Replace linear weights with neural networks
2. **Implement PSO**: Particle swarm optimization for continuous spaces
3. **Load Real PDB Data**: Train on actual protein structures
4. **Advanced Functions**: Secondary structure prediction, loop modeling
5. **Visualization**: 3D structure visualization and analysis

## Future Enhancements

- **Multi-objective Evolution**: Optimize multiple fitness criteria simultaneously
- **Parallel Processing**: GPU acceleration and distributed computing
- **Web Interface**: Browser-based parameter tuning and visualization
- **Advanced Evaluation**: More sophisticated fitness functions

## License

This project is released under the MIT License.

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Submit a pull request

---

*Built with modern C++17 for performance and maintainability.*
