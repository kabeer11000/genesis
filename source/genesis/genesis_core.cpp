#include "genesis_core.h"
#include "genesis_config.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <filesystem>
#include <cstring>
#include <random> // Required for random number generation
#include <chrono>   // Required for seeding the random number generator

// Amount of functions to generate
#define N 50
// Amino acid properties (Kyte-Doolittle hydrophobicity scale)
static std::unordered_map<char, std::vector<double>> aa_properties = {
    {'A', {1.8, 67.0, 0.0, 71.08}},     // Ala: hydrophobicity, volume, charge, mass
    {'R', {-4.5, 148.0, 1.0, 156.19}},  // Arg
    {'N', {-3.5, 96.0, 0.0, 114.10}},   // Asn
    {'D', {-3.5, 91.0, -1.0, 115.09}},  // Asp
    {'C', {2.5, 86.0, 0.0, 103.14}},    // Cys
    {'Q', {-3.5, 114.0, 0.0, 128.13}},  // Gln
    {'E', {-3.5, 109.0, -1.0, 129.12}}, // Glu
    {'G', {-0.4, 48.0, 0.0, 57.05}},    // Gly
    {'H', {-3.2, 118.0, 0.5, 137.14}},  // His
    {'I', {4.5, 124.0, 0.0, 113.16}},   // Ile
    {'L', {3.8, 124.0, 0.0, 113.16}},   // Leu
    {'K', {-3.9, 135.0, 1.0, 128.17}},  // Lys
    {'M', {1.9, 124.0, 0.0, 131.20}},   // Met
    {'F', {2.8, 135.0, 0.0, 147.18}},   // Phe
    {'P', {-1.6, 90.0, 0.0, 97.12}},    // Pro
    {'S', {-0.8, 73.0, 0.0, 87.08}},    // Ser
    {'T', {-0.7, 93.0, 0.0, 101.11}},   // Thr
    {'W', {-0.9, 163.0, 0.0, 186.21}},  // Trp
    {'Y', {-1.3, 141.0, 0.0, 163.18}},  // Tyr
    {'V', {4.2, 105.0, 0.0, 99.13}}     // Val
};

// AminoAcid implementation
AminoAcid::AminoAcid(char c) : code(c)
{
    if (aa_properties.find(c) != aa_properties.end())
    {
        auto props = aa_properties[c];
        hydrophobicity = props[0];
        volume = props[1];
        charge = props[2];
        mass = props[3];
    }
    else
    {
        // Default values for unknown amino acids
        hydrophobicity = 0.0;
        volume = 100.0;
        charge = 0.0;
        mass = 100.0;
    }
}

// Protein implementation
Protein::Protein(const std::string &seq) : sequence(seq)
{
    amino_acids.reserve(seq.length());
    for (char c : seq)
    {
        amino_acids.emplace_back(c);
    }

    // Initialize coordinate vectors
    backbone_coords.resize(seq.length());
    all_atom_coords.resize(seq.length() * 4); // Approximate 4 heavy atoms per residue
    phi_angles.resize(seq.length(), 0.0);
    psi_angles.resize(seq.length(), 0.0);
    local_features.resize(seq.length() * 10, 0.0); // Extra space for function outputs
}

Protein Protein::fromPDB(const std::string &filename)
{
    std::ifstream file(filename);
    std::string line;
    std::vector<std::string> sequence_chars;
    std::vector<Coord3D> ca_coords;

    while (std::getline(file, line))
    {
        if (line.substr(0, 4) == "ATOM" && line.substr(12, 4) == " CA ")
        {
            // Extract CA coordinates
            double x = std::stod(line.substr(30, 8));
            double y = std::stod(line.substr(38, 8));
            double z = std::stod(line.substr(46, 8));
            ca_coords.emplace_back(x, y, z);

            // Extract amino acid
            std::string res_name = line.substr(17, 3);
            char aa_code = 'A'; // Default

            // Convert 3-letter to 1-letter code (simplified)
            if (res_name == "ALA")
                aa_code = 'A';
            else if (res_name == "ARG")
                aa_code = 'R';
            else if (res_name == "ASN")
                aa_code = 'N';
            else if (res_name == "ASP")
                aa_code = 'D';
            else if (res_name == "CYS")
                aa_code = 'C';
            else if (res_name == "GLN")
                aa_code = 'Q';
            else if (res_name == "GLU")
                aa_code = 'E';
            else if (res_name == "GLY")
                aa_code = 'G';
            else if (res_name == "HIS")
                aa_code = 'H';
            else if (res_name == "ILE")
                aa_code = 'I';
            else if (res_name == "LEU")
                aa_code = 'L';
            else if (res_name == "LYS")
                aa_code = 'K';
            else if (res_name == "MET")
                aa_code = 'M';
            else if (res_name == "PHE")
                aa_code = 'F';
            else if (res_name == "PRO")
                aa_code = 'P';
            else if (res_name == "SER")
                aa_code = 'S';
            else if (res_name == "THR")
                aa_code = 'T';
            else if (res_name == "TRP")
                aa_code = 'W';
            else if (res_name == "TYR")
                aa_code = 'Y';
            else if (res_name == "VAL")
                aa_code = 'V';

            sequence_chars.push_back(std::string(1, aa_code));
        }
    }

    // Build sequence string
    std::string sequence;
    for (const auto &aa : sequence_chars)
    {
        sequence += aa;
    }

    Protein protein(sequence);
    protein.backbone_coords = ca_coords;

    return protein;
}

double Protein::calculateRMSD(const Protein &other) const
{
    if (backbone_coords.size() != other.backbone_coords.size())
    {
        return 1000.0; // Large penalty for size mismatch
    }

    double sum_sq_dist = 0.0;
    for (size_t i = 0; i < backbone_coords.size(); ++i)
    {
        double dist = backbone_coords[i].distance(other.backbone_coords[i]);
        sum_sq_dist += dist * dist;
    }

    return sqrt(sum_sq_dist / backbone_coords.size());
}

const AminoAcid &Protein::getAminoAcid(size_t pos) const
{
    return amino_acids[pos];
}

void Protein::updateCoordinates(const std::vector<Coord3D> &coords)
{
    backbone_coords = coords;
}

// Function implementations
void AngleFunction::apply(Protein &protein, const std::vector<double> &weights)
{
    double phi_weight = weights[0];
    double psi_weight = weights[1];

    // Calculate dihedral angles and update coordinates based on them
    for (size_t i = 1; i < protein.backbone_coords.size() - 1; ++i)
    {
        // Simplified angle calculation
        Coord3D &prev = protein.backbone_coords[i - 1];
        Coord3D &curr = protein.backbone_coords[i];
        Coord3D &next = protein.backbone_coords[i + 1];

        // Calculate phi angle (simplified)
        double phi = atan2(curr.y - prev.y, curr.x - prev.x);
        double psi = atan2(next.y - curr.y, next.x - curr.x);

        protein.phi_angles[i] = phi;
        protein.psi_angles[i] = psi;

        // Adjust coordinates based on weighted angles
        curr.x += phi_weight * cos(phi) * 0.1;
        curr.y += phi_weight * sin(phi) * 0.1;
        curr.z += psi_weight * sin(psi) * 0.1;
    }
}

void HydrophobicityFunction::apply(Protein &protein, const std::vector<double> &weights)
{
    double x_weight = weights[0];
    double y_weight = weights[1];
    double z_weight = weights[2];

    // Adjust coordinates based on hydrophobicity
    for (size_t i = 0; i < protein.backbone_coords.size(); ++i)
    {
        double hydro = protein.getAminoAcid(i).hydrophobicity;

        protein.backbone_coords[i].x += x_weight * hydro * 0.05;
        protein.backbone_coords[i].y += y_weight * hydro * 0.05;
        protein.backbone_coords[i].z += z_weight * hydro * 0.05;
    }
}

void CombinationFunction::apply(Protein &protein, const std::vector<double> &weights)
{
    // Combine previous function results stored in local_features
    for (size_t i = 0; i < protein.backbone_coords.size(); ++i)
    {
        Coord3D adjustment(0, 0, 0);

        for (size_t j = 0; j < input_function_indices.size() && j < weights.size(); ++j)
        {
            size_t feature_idx = input_function_indices[j] * protein.length() + i;
            if (feature_idx < protein.local_features.size())
            {
                double feature_value = protein.local_features[feature_idx];
                adjustment.x += weights[j] * feature_value * 0.01;
                adjustment.y += weights[j] * feature_value * 0.01;
                adjustment.z += weights[j] * feature_value * 0.01;
            }
        }

        protein.backbone_coords[i] = protein.backbone_coords[i] + adjustment;
    }
}

void DistanceConstraintFunction::apply(Protein &protein, const std::vector<double> &weights)
{
    double min_dist = weights[0];
    double max_dist = weights[1];
    double strength = weights[2];
    double range = std::max(1.0, weights[3]);

    // Apply distance constraints between nearby residues
    for (size_t i = 0; i < protein.backbone_coords.size(); ++i)
    {
        for (size_t j = i + 1; j < std::min(protein.backbone_coords.size(),
                                            i + static_cast<size_t>(range));
             ++j)
        {
            double dist = protein.backbone_coords[i].distance(protein.backbone_coords[j]);

            if (dist < min_dist)
            {
                // Push apart
                Coord3D direction = protein.backbone_coords[j] +
                                    (protein.backbone_coords[i] * -1.0);
                double factor = strength * (min_dist - dist) / (dist + 1e-6);
                protein.backbone_coords[j] = protein.backbone_coords[j] +
                                             (direction * (factor * 0.5));
            }
            else if (dist > max_dist)
            {
                // Pull together
                Coord3D direction = protein.backbone_coords[i] +
                                    (protein.backbone_coords[j] * -1.0);
                double factor = strength * (dist - max_dist) / (dist + 1e-6);
                protein.backbone_coords[j] = protein.backbone_coords[j] +
                                             (direction * (factor * 0.5));
            }
        }
    }
}

// Chromosome implementation
Chromosome::Chromosome(size_t weight_count) : fitness(0.0)
{
    weights.resize(weight_count);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, 1.0);

    for (double &w : weights)
    {
        w = dist(gen);
    }
}

Chromosome Chromosome::crossover(const Chromosome &other) const
{
    std::vector<double> child_weights(weights.size());
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (size_t i = 0; i < weights.size(); ++i)
    {
        if (dist(gen) < 0.5)
        {
            child_weights[i] = weights[i];
        }
        else
        {
            child_weights[i] = other.weights[i];
        }
    }

    return Chromosome(child_weights);
}

void Chromosome::mutate(double mutation_rate, double mutation_strength)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> rate_dist(0.0, 1.0);
    std::normal_distribution<double> strength_dist(0.0, mutation_strength);

    for (double &w : weights)
    {
        if (rate_dist(gen) < mutation_rate)
        {
            w += strength_dist(gen);
        }
    }
}

// FunctionExecutor implementation
void FunctionExecutor::addFunction(std::unique_ptr<Function> func)
{
    if (!functions.empty())
    {
        weight_offsets.push_back(weight_offsets.back() + functions.back()->getWeightCount());
    }
    else
    {
        weight_offsets.push_back(0);
    }
    functions.push_back(std::move(func));
}

Protein FunctionExecutor::execute(const Protein &input, const Chromosome &chromosome)
{
    Protein result = input; // Copy

    for (size_t i = 0; i < functions.size(); ++i)
    {
        size_t start_idx = weight_offsets[i];
        size_t weight_count = functions[i]->getWeightCount();

        std::vector<double> func_weights(
            chromosome.weights.begin() + start_idx,
            chromosome.weights.begin() + start_idx + weight_count);

        functions[i]->apply(result, func_weights);
    }

    return result;
}

size_t FunctionExecutor::getTotalWeightCount() const
{
    size_t total = 0;
    for (const auto &func : functions)
    {
        total += func->getWeightCount();
    }
    return total;
}

FunctionExecutor FunctionExecutor::createStandardSet()
{
    FunctionExecutor executor;

    // Add predefined functions
    executor.addFunction(std::make_unique<AngleFunction>());
    executor.addFunction(std::make_unique<HydrophobicityFunction>());
    executor.addFunction(std::make_unique<DistanceConstraintFunction>());

    // Add combination functions that use previous functions
    executor.addFunction(std::make_unique<CombinationFunction>(std::vector<size_t>{0, 1}));
    executor.addFunction(std::make_unique<CombinationFunction>(std::vector<size_t>{1, 2}));
    executor.addFunction(std::make_unique<CombinationFunction>(std::vector<size_t>{0, 2, 3}));

    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());

    for (int i = 6; i < N; ++i)
    {
        // Determine the number of dependencies you want for this iteration.
        // Let's say between 1 and 4 dependencies for example.
        std::uniform_int_distribution<size_t> num_deps_dist(1, 5);
        size_t num_dependencies = num_deps_dist(rng);

        std::vector<size_t> dependencies;
        dependencies.reserve(num_dependencies); // Pre-allocate memory

        // Generate random dependencies, each smaller than 'i'
        // Ensure 'i' is at least 1 for the uniform_int_distribution to be valid.
        // If i is 0 or 1, max value will be 0, which is fine, but practically you'll likely have higher 'i'.
        std::uniform_int_distribution<size_t> dep_value_dist(0, i - 1);

        for (size_t k = 0; k < num_dependencies; ++k)
        {
            dependencies.push_back(dep_value_dist(rng));
        }

        // Optional: If you want unique dependencies, you can sort and unique them,
        // though this might reduce the number of dependencies if duplicates are generated.
        std::sort(dependencies.begin(), dependencies.end());
        dependencies.erase(std::unique(dependencies.begin(), dependencies.end()), dependencies.end());

        executor.addFunction(std::make_unique<CombinationFunction>(dependencies));
    }

    return executor;
}

// FitnessEvaluator implementation
void FitnessEvaluator::addReference(const Protein &protein)
{
    reference_proteins.push_back(protein);
}

void FitnessEvaluator::loadReferencesFromPDBDirectory(const std::string &directory)
{
    try
    {
        for (const auto &entry : std::filesystem::directory_iterator(directory))
        {
            if (entry.path().extension() == ".pdb")
            {
                Protein ref = Protein::fromPDB(entry.path().string());
                reference_proteins.push_back(ref);
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << "Error loading PDB files: " << e.what() << std::endl;
    }
}

double FitnessEvaluator::evaluate(const Protein &predicted, const std::string &target_sequence)
{
    if (predicted.sequence != target_sequence)
    {
        return 0.0; // Sequence mismatch
    }

    double best_fitness = 0.0;

    // Find the best matching reference protein
    for (const auto &ref : reference_proteins)
    {
        if (ref.sequence == target_sequence)
        {
            double rmsd = predicted.calculateRMSD(ref);
            double fitness = 1.0 / (1.0 + rmsd); // Convert RMSD to fitness
            best_fitness = std::max(best_fitness, fitness);
        }
    }

    // If no exact sequence match, use structural features
    if (best_fitness == 0.0)
    {
        best_fitness = evaluateStructuralFeatures(predicted);
    }

    return best_fitness;
}

double FitnessEvaluator::evaluateStructuralFeatures(const Protein &predicted)
{
    double fitness = 0.0;

    // Evaluate based on reasonable structural constraints
    // 1. Bond length consistency
    double bond_score = 0.0;
    for (size_t i = 1; i < predicted.backbone_coords.size(); ++i)
    {
        double dist = predicted.backbone_coords[i - 1].distance(predicted.backbone_coords[i]);
        // Ideal CA-CA distance is ~3.8 Angstroms
        double ideal_dist = 3.8;
        bond_score += 1.0 / (1.0 + abs(dist - ideal_dist));
    }
    bond_score /= (predicted.backbone_coords.size() - 1);

    // 2. Compactness (avoiding extended conformations)
    double compactness = 0.0;
    Coord3D center(0, 0, 0);
    for (const auto &coord : predicted.backbone_coords)
    {
        center = center + coord;
    }
    center = center * (1.0 / predicted.backbone_coords.size());

    double avg_distance = 0.0;
    for (const auto &coord : predicted.backbone_coords)
    {
        avg_distance += coord.distance(center);
    }
    avg_distance /= predicted.backbone_coords.size();

    compactness = 1.0 / (1.0 + avg_distance / 10.0); // Normalize

    fitness = 0.7 * bond_score + 0.3 * compactness;

    return fitness;
}

// EvolutionaryAlgorithm implementation
EvolutionaryAlgorithm::EvolutionaryAlgorithm(size_t pop_size, double mut_rate,
                                             double mut_strength, double cross_rate)
    : population_size(pop_size), mutation_rate(mut_rate),
      mutation_strength(mut_strength), crossover_rate(cross_rate),
      executor(FunctionExecutor::createStandardSet())
{
    std::random_device rd;
    rng.seed(rd());
}

void EvolutionaryAlgorithm::initialize(size_t weight_count)
{
    population.clear();
    population.reserve(population_size);

    for (size_t i = 0; i < population_size; ++i)
    {
        population.emplace_back(weight_count);
    }
}

void EvolutionaryAlgorithm::evolve(size_t generations, const std::vector<std::string> &target_sequences)
{
    for (size_t gen = 0; gen < generations; ++gen)
    {
        // Evaluate fitness
        for (auto &chromosome : population)
        {
            double total_fitness = 0.0;

            for (const auto &target_seq : target_sequences)
            {
                Protein input(target_seq);
                Protein predicted = executor.execute(input, chromosome);
                double fitness = evaluator.evaluate(predicted, target_seq);
                total_fitness += fitness;
            }

            chromosome.fitness = total_fitness / target_sequences.size();
        }

        // Sort by fitness
        std::sort(population.begin(), population.end());

        // Print progress
        if (gen % 10 == 0)
        {
            std::cout << "Generation " << gen << ", Best fitness: "
                      << population[0].fitness << std::endl;
        }

        // Create next generation
        std::vector<Chromosome> next_generation;

        // Keep top 10% (elitism)
        size_t elite_count = population_size / 10;
        for (size_t i = 0; i < elite_count; ++i)
        {
            next_generation.push_back(population[i]);
        }

        // Generate rest through crossover and mutation
        std::uniform_real_distribution<double> dist(0.0, 1.0);
        std::uniform_int_distribution<size_t> parent_dist(0, population_size / 2);

        while (next_generation.size() < population_size)
        {
            size_t parent1_idx = parent_dist(rng);
            size_t parent2_idx = parent_dist(rng);

            Chromosome child = population[parent1_idx].crossover(population[parent2_idx]);
            child.mutate(mutation_rate, mutation_strength);
            next_generation.push_back(child);
        }

        population = next_generation;
    }
}

Chromosome EvolutionaryAlgorithm::getBest() const
{
    return *std::max_element(population.begin(), population.end(),
                             [](const Chromosome &a, const Chromosome &b)
                             {
                                 return a.fitness < b.fitness;
                             });
}

std::vector<double> EvolutionaryAlgorithm::getPopulationFitness() const
{
    std::vector<double> fitness_values;
    for (const auto &chromosome : population)
    {
        fitness_values.push_back(chromosome.fitness);
    }
    return fitness_values;
}

void EvolutionaryAlgorithm::addReferenceProtein(const Protein &protein)
{
    evaluator.addReference(protein);
}

void EvolutionaryAlgorithm::loadReferencePDBs(const std::string &directory)
{
    evaluator.loadReferencesFromPDBDirectory(directory);
}

// ConfigurableExperiment implementation
void ConfigurableExperiment::setupExperiment()
{
    // Create evolutionary algorithm with config parameters
    ea = std::make_unique<EvolutionaryAlgorithm>(
        config.population_size,
        config.mutation_rate,
        config.mutation_strength,
        config.crossover_rate);

    // Setup function executor based on config
    executor = std::make_unique<FunctionExecutor>();

    if (config.use_angle_function)
    {
        executor->addFunction(std::make_unique<AngleFunction>());
    }

    if (config.use_hydrophobicity_function)
    {
        executor->addFunction(std::make_unique<HydrophobicityFunction>());
    }

    if (config.use_distance_constraint_function)
    {
        executor->addFunction(std::make_unique<DistanceConstraintFunction>());
    }

    if (config.use_combination_functions)
    {
        // Add combination functions that use previous functions
        executor->addFunction(std::make_unique<CombinationFunction>(std::vector<size_t>{0, 1}));
        executor->addFunction(std::make_unique<CombinationFunction>(std::vector<size_t>{1, 2}));
    }
}

void ConfigurableExperiment::run()
{
    if (config.verbose)
    {
        std::cout << "=== Configurable Genesis Experiment ===" << std::endl;
        printConfig();
    }

    // Initialize population
    size_t weight_count = executor->getTotalWeightCount();
    ea->initialize(weight_count);

    if (config.verbose)
    {
        std::cout << "Chromosome size: " << weight_count << " weights" << std::endl;
        std::cout << "Starting evolution..." << std::endl;
    }

    // Run evolution
    auto start_time = std::chrono::high_resolution_clock::now();
    ea->evolve(config.generations, config.target_sequences);
    auto end_time = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    if (config.verbose)
    {
        printResults(duration.count());
    }

    if (config.save_results)
    {
        saveResults();
    }
}

void ConfigurableExperiment::printConfig()
{
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Population size: " << config.population_size << std::endl;
    std::cout << "  Generations: " << config.generations << std::endl;
    std::cout << "  Mutation rate: " << config.mutation_rate << std::endl;
    std::cout << "  Mutation strength: " << config.mutation_strength << std::endl;
    std::cout << "  Crossover rate: " << config.crossover_rate << std::endl;
    std::cout << "  Target sequences: " << config.target_sequences.size() << std::endl;

    for (size_t i = 0; i < config.target_sequences.size(); ++i)
    {
        std::cout << "    " << i + 1 << ": " << config.target_sequences[i]
                  << " (" << config.target_sequences[i].length() << " residues)" << std::endl;
    }
    std::cout << std::endl;
}

void ConfigurableExperiment::printResults(long long duration_ms)
{
    Chromosome best = ea->getBest();
    std::vector<double> fitness_values = ea->getPopulationFitness();

    double sum_fitness = 0.0;
    double max_fitness = 0.0;
    double min_fitness = 1000.0;

    for (double f : fitness_values)
    {
        sum_fitness += f;
        max_fitness = std::max(max_fitness, f);
        min_fitness = std::min(min_fitness, f);
    }
    double avg_fitness = sum_fitness / fitness_values.size();

    std::cout << "=== Results ===" << std::endl;
    std::cout << "Best fitness: " << std::fixed << std::setprecision(6) << best.fitness << std::endl;
    std::cout << "Average fitness: " << avg_fitness << std::endl;
    std::cout << "Min fitness: " << min_fitness << std::endl;
    std::cout << "Max fitness: " << max_fitness << std::endl;
    std::cout << "Evolution time: " << duration_ms << " ms" << std::endl;
    std::cout << "Time per generation: " << duration_ms / config.generations << " ms" << std::endl;
    std::cout << std::endl;
}

void ConfigurableExperiment::saveResults()
{
    // TODO: Implement result saving to files
    std::cout << "Result saving not yet implemented" << std::endl;
}

Chromosome ConfigurableExperiment::getBestChromosome()
{
    return ea->getBest();
}

std::vector<double> ConfigurableExperiment::getPopulationFitness()
{
    return ea->getPopulationFitness();
}

// BatchExperimentRunner implementation
void BatchExperimentRunner::runMutationRateSweep()
{
    std::cout << "=== Mutation Rate Sweep ===" << std::endl;
    std::vector<double> mutation_rates = {0.01, 0.05, 0.1, 0.15, 0.2, 0.3};

    std::cout << std::setw(12) << "MutationRate" << std::setw(12) << "BestFitness"
              << std::setw(12) << "AvgFitness" << std::setw(10) << "Time(ms)" << std::endl;
    std::cout << std::string(46, '-') << std::endl;

    for (double rate : mutation_rates)
    {
        ExperimentConfig config = ExperimentConfigs::getParameterSweep();
        config.mutation_rate = rate;
        config.verbose = false;

        ConfigurableExperiment exp(config);

        auto start_time = std::chrono::high_resolution_clock::now();
        exp.run();
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        Chromosome best = exp.getBestChromosome();
        std::vector<double> fitness = exp.getPopulationFitness();
        double avg_fitness = 0.0;
        for (double f : fitness)
            avg_fitness += f;
        avg_fitness /= fitness.size();

        std::cout << std::setw(12) << std::fixed << std::setprecision(3) << rate
                  << std::setw(12) << std::setprecision(6) << best.fitness
                  << std::setw(12) << avg_fitness
                  << std::setw(10) << duration.count() << std::endl;
    }
    std::cout << std::endl;
}

void BatchExperimentRunner::runPopulationSizeSweep()
{
    std::cout << "=== Population Size Sweep ===" << std::endl;
    std::vector<size_t> population_sizes = {20, 50, 100, 150, 200};

    std::cout << std::setw(12) << "PopSize" << std::setw(12) << "BestFitness"
              << std::setw(12) << "AvgFitness" << std::setw(10) << "Time(ms)" << std::endl;
    std::cout << std::string(46, '-') << std::endl;

    for (size_t pop_size : population_sizes)
    {
        ExperimentConfig config = ExperimentConfigs::getParameterSweep();
        config.population_size = pop_size;
        config.verbose = false;

        ConfigurableExperiment exp(config);

        auto start_time = std::chrono::high_resolution_clock::now();
        exp.run();
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        Chromosome best = exp.getBestChromosome();
        std::vector<double> fitness = exp.getPopulationFitness();
        double avg_fitness = 0.0;
        for (double f : fitness)
            avg_fitness += f;
        avg_fitness /= fitness.size();

        std::cout << std::setw(12) << pop_size
                  << std::setw(12) << std::fixed << std::setprecision(6) << best.fitness
                  << std::setw(12) << avg_fitness
                  << std::setw(10) << duration.count() << std::endl;
    }
    std::cout << std::endl;
}

void BatchExperimentRunner::runFunctionConfigurationTest()
{
    std::cout << "=== Function Configuration Test ===" << std::endl;

    struct ConfigTest
    {
        std::string name;
        bool angle, hydro, distance, combination;
    };

    std::vector<ConfigTest> configs = {
        {"Angles Only", true, false, false, false},
        {"Hydro Only", false, true, false, false},
        {"Distance Only", false, false, true, false},
        {"Angle+Hydro", true, true, false, false},
        {"All Basic", true, true, true, false},
        {"All+Combination", true, true, true, true}};

    std::cout << std::setw(15) << "Configuration" << std::setw(12) << "BestFitness"
              << std::setw(10) << "Weights" << std::setw(10) << "Time(ms)" << std::endl;
    std::cout << std::string(47, '-') << std::endl;

    for (const auto &test : configs)
    {
        ExperimentConfig config = ExperimentConfigs::getMiniTest();
        config.use_angle_function = test.angle;
        config.use_hydrophobicity_function = test.hydro;
        config.use_distance_constraint_function = test.distance;
        config.use_combination_functions = test.combination;
        config.verbose = false;

        ConfigurableExperiment exp(config);

        auto start_time = std::chrono::high_resolution_clock::now();
        exp.run();
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

        Chromosome best = exp.getBestChromosome();

        std::cout << std::setw(15) << test.name
                  << std::setw(12) << std::fixed << std::setprecision(6) << best.fitness
                  << std::setw(10) << best.weights.size()
                  << std::setw(10) << duration.count() << std::endl;
    }
    std::cout << std::endl;
}
