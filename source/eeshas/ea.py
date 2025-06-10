import numpy as np
from features import FEATURE_FUNCTIONS
from utils import rmsd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class EvolutionaryAlgorithm:
    def __init__(self, proteins, population_size=100, generations=50, mutation_rate=0.1):
        self.proteins = proteins
        self.population_size = population_size
        self.generations = generations
        self.mutation_rate = mutation_rate
        self.n_features = len(FEATURE_FUNCTIONS) * 10
        
        # Initialize PCA for dimensionality reduction
        self.pca = PCA(n_components=3)
        self.scaler = StandardScaler()
        self._prepare_feature_space()
    
    def _prepare_feature_space(self):
        """Prepare feature extraction and transformation"""
        all_features = []
        for protein in self.proteins:
            features = self._extract_features(protein)
            all_features.append(features)
        
        # Fit scaler and PCA
        all_features = np.array(all_features)
        scaled_features = self.scaler.fit_transform(all_features)
        self.pca.fit(scaled_features)
    
    def _extract_features(self, protein):
        """Extract all features from a protein"""
        features = []
        for func in FEATURE_FUNCTIONS:
            features.extend(func(protein))
        return np.array(features)
    
    def generate_individual(self):
        """Generate a random individual (transformation matrix)"""
        # Individual represents a linear transformation from feature space to 3D space
        return np.random.randn(self.n_features, 3) * 0.1
    
    def features_to_coords(self, protein, transform_matrix):
        """Convert protein features to 3D coordinates using transformation matrix"""
        n_residues = len(protein['coords'])
        
        # Extract features
        features = self._extract_features(protein)
        
        # Scale features
        features_scaled = self.scaler.transform([features])[0]
        
        # Apply transformation to get base 3D representation
        base_coords = np.dot(features_scaled, transform_matrix)
        
        # Generate coordinates for each residue
        # Use a simple linear interpolation along the protein chain
        coords = np.zeros((n_residues, 3))
        for i in range(n_residues):
            # Create position-dependent variation
            t = i / max(1, n_residues - 1)
            # Interpolate and add position-dependent transformation
            coords[i] = base_coords * (1 + 0.1 * np.sin(2 * np.pi * t * np.array([1, 2, 3])))
            coords[i] += np.array([
                10 * np.cos(0.1 * i),
                10 * np.sin(0.1 * i),
                0.5 * i
            ])
        
        return coords
    
    def evaluate(self, individual):
        """Evaluate fitness of an individual"""
        total_rmsd = 0
        count = 0
        
        for protein in self.proteins:
            try:
                # Generate predicted coordinates
                pred_coords = self.features_to_coords(protein, individual)
                true_coords = protein['coords']
                print("true_coords: ", true_coords, "pred_coords: ", pred_coords)
                # Calculate RMSD
                if pred_coords.shape == true_coords.shape:
                    rmsd_value = rmsd(pred_coords, true_coords)
                    total_rmsd += rmsd_value
                    count += 1
            except Exception as e:
                # Skip problematic proteins
                continue
        
        # Return average RMSD (lower is better)
        return total_rmsd / max(1, count)
    
    def crossover(self, parent1, parent2):
        """Perform crossover between two parents"""
        # Uniform crossover at the gene level
        mask = np.random.rand(*parent1.shape) < 0.5
        child = np.where(mask, parent1, parent2)
        return child
    
    def mutate(self, individual):
        """Mutate an individual"""
        mutation_mask = np.random.rand(*individual.shape) < self.mutation_rate
        mutations = np.random.randn(*individual.shape) * 0.1
        individual[mutation_mask] += mutations[mutation_mask]
        return individual
    
    def run(self):
        """Run the evolutionary algorithm"""
        # Initialize population
        population = [self.generate_individual() for _ in range(self.population_size)]
        
        best_fitness_history = []
        
        for gen in range(self.generations):
            # Evaluate fitness for all individuals
            fitness_scores = []
            for ind in population:
                fitness = self.evaluate(ind)
                fitness_scores.append(fitness)
            
            # Sort by fitness (lower RMSD is better)
            sorted_indices = np.argsort(fitness_scores)
            sorted_population = [population[i] for i in sorted_indices]
            sorted_fitness = [fitness_scores[i] for i in sorted_indices]
            
            # Print progress
            best_fitness = sorted_fitness[0]
            avg_fitness = np.mean(sorted_fitness)
            print(f"Generation {gen + 1}/{self.generations}")
            print(f"  Best RMSD: {best_fitness:.3f}")
            print(f"  Avg RMSD:  {avg_fitness:.3f}")
            best_fitness_history.append(best_fitness)
            
            # Selection and reproduction
            new_population = []
            
            # Elitism - keep best individual
            new_population.append(sorted_population[0].copy())
            
            # Generate rest of population
            while len(new_population) < self.population_size:
                # Tournament selection
                tournament_size = 3
                tournament_indices = np.random.choice(
                    len(sorted_population) // 2,  # Select from better half
                    size=tournament_size,
                    replace=False
                )
                parent1_idx = tournament_indices[np.argmin([sorted_fitness[i] for i in tournament_indices])]
                
                tournament_indices = np.random.choice(
                    len(sorted_population) // 2,
                    size=tournament_size,
                    replace=False
                )
                parent2_idx = tournament_indices[np.argmin([sorted_fitness[i] for i in tournament_indices])]
                
                # Crossover
                child = self.crossover(sorted_population[parent1_idx], sorted_population[parent2_idx])
                
                # Mutation
                child = self.mutate(child)
                
                new_population.append(child)
            
            population = new_population
        
        # Return best individual
        final_fitness = [self.evaluate(ind) for ind in population]
        best_idx = np.argmin(final_fitness)
        
        print("\nEvolution complete!")
        print(f"Final best RMSD: {final_fitness[best_idx]:.3f}")
        
        return population[best_idx], best_fitness_history
