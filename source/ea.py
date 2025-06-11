# import numpy as np
# from features import FEATURE_FUNCTIONS
# from utils import rmsd
# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler
# import matplotlib.pyplot as plt
# import pandas as pd
# import time
# import logging
# from datetime import datetime


# class EvolutionaryAlgorithm:
#     def __init__(self, proteins, population_size=100, generations=50, mutation_rate=0.1, 
#                  selection_method='tournament', crossover_method='uniform', 
#                  mutation_method='gaussian', log_level='INFO'):
#         self.proteins = proteins
#         self.population_size = population_size
#         self.generations = generations
#         self.mutation_rate = mutation_rate
#         self.selection_method = selection_method
#         self.crossover_method = crossover_method
#         self.mutation_method = mutation_method
#         self.n_features = len(FEATURE_FUNCTIONS) * 10
        
#         # Setup logging first
#         self.setup_logging(log_level)
        
#         # Initialize PCA for dimensionality reduction
#         self.pca = PCA(n_components=3)
#         self.scaler = StandardScaler()
#         self._prepare_feature_space()
#         self.stats = {
#             'fitness_history': [],
#             'avg_fitness_history': [],
#             'diversity_history': [],
#             'selection_pressure': [],
#             'convergence_gen': None,
#             'execution_time': 0
#         }
    
#     def setup_logging(self, level):
#         """Setup logging configuration"""
#         logging.basicConfig(
#             level=getattr(logging, level),
#             format='%(asctime)s - %(levelname)s - %(message)s',
#             handlers=[
#                 logging.FileHandler(f'ea_log_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'),
#                 logging.StreamHandler()
#             ]
#         )
#         self.logger = logging.getLogger(__name__)
    
#     def _prepare_feature_space(self):
#         """Prepare feature extraction and transformation"""
#         all_features = []
#         for protein in self.proteins:
#             features = self._extract_features(protein)
#             all_features.append(features)
        
#         # Fit scaler and PCA
#         all_features = np.array(all_features)
#         scaled_features = self.scaler.fit_transform(all_features)
#         self.pca.fit(scaled_features)
#         self.logger.info(f"Feature space prepared: {len(all_features)} proteins, {self.n_features} features")
    
#     def _extract_features(self, protein):
#         """Extract all features from a protein"""
#         features = []
#         for func in FEATURE_FUNCTIONS:
#             features.extend(func(protein))
#         return np.array(features)
    
#     def generate_individual(self):
#         """Generate a random individual (transformation matrix)"""
#         return np.random.randn(self.n_features, 3) * 0.1
    
#     def features_to_coords(self, protein, transform_matrix):
#         """Convert protein features to 3D coordinates using transformation matrix"""
#         n_residues = len(protein['coords'])
        
#         # Extract features
#         features = self._extract_features(protein)
        
#         # Scale features
#         features_scaled = self.scaler.transform([features])[0]
        
#         # Apply transformation to get base 3D representation
#         base_coords = np.dot(features_scaled, transform_matrix)
        
#         # Generate coordinates for each residue
#         coords = np.zeros((n_residues, 3))
#         for i in range(n_residues):
#             t = i / max(1, n_residues - 1)
#             coords[i] = base_coords * (1 + 0.1 * np.sin(2 * np.pi * t * np.array([1, 2, 3])))
#             coords[i] += np.array([
#                 10 * np.cos(0.1 * i),
#                 10 * np.sin(0.1 * i),
#                 0.5 * i
#             ])
        
#         return coords
    
#     def evaluate(self, individual):
#         """Evaluate fitness of an individual"""
#         total_rmsd = 0
#         count = 0
        
#         for protein in self.proteins:
#             try:
#                 pred_coords = self.features_to_coords(protein, individual)
#                 true_coords = protein['coords']
                
#                 if pred_coords.shape == true_coords.shape:
#                     rmsd_value = rmsd(pred_coords, true_coords)
#                     total_rmsd += rmsd_value
#                     count += 1
#             except Exception as e:
#                 self.logger.warning(f"Error evaluating protein: {e}")
#                 continue
        
#         return total_rmsd / max(1, count)
    
#     def tournament_selection(self, population, fitness_scores, tournament_size=3):
#         """Tournament selection"""
#         tournament_indices = np.random.choice(len(population), size=tournament_size, replace=False)
#         winner_idx = tournament_indices[np.argmin([fitness_scores[i] for i in tournament_indices])]
#         return population[winner_idx]
    
#     def binary_tournament_selection(self, population, fitness_scores):
#         """Binary tournament selection"""
#         return self.tournament_selection(population, fitness_scores, tournament_size=2)
    
#     def roulette_wheel_selection(self, population, fitness_scores):
#         """Roulette wheel selection (fitness proportionate)"""
#         # Convert to maximization problem (lower RMSD is better)
#         max_fitness = max(fitness_scores)
#         weights = [max_fitness - f + 1e-10 for f in fitness_scores]
#         weights = np.array(weights) / sum(weights)
#         selected_idx = np.random.choice(len(population), p=weights)
#         return population[selected_idx]
    
#     def rank_selection(self, population, fitness_scores):
#         """Rank-based selection"""
#         sorted_indices = np.argsort(fitness_scores)
#         ranks = np.zeros(len(population))
#         ranks[sorted_indices] = np.arange(len(population), 0, -1)
#         weights = ranks / sum(ranks)
#         selected_idx = np.random.choice(len(population), p=weights)
#         return population[selected_idx]
    
#     def select_parent(self, population, fitness_scores):
#         """Select parent based on selection method"""
#         if self.selection_method == 'tournament':
#             return self.tournament_selection(population, fitness_scores)
#         elif self.selection_method == 'binary_tournament':
#             return self.binary_tournament_selection(population, fitness_scores)
#         elif self.selection_method == 'roulette':
#             return self.roulette_wheel_selection(population, fitness_scores)
#         elif self.selection_method == 'rank':
#             return self.rank_selection(population, fitness_scores)
#         else:
#             raise ValueError(f"Unknown selection method: {self.selection_method}")
    
#     def uniform_crossover(self, parent1, parent2):
#         """Uniform crossover"""
#         mask = np.random.rand(*parent1.shape) < 0.5
#         return np.where(mask, parent1, parent2)
    
#     def single_point_crossover(self, parent1, parent2):
#         """Single point crossover"""
#         point = np.random.randint(1, parent1.shape[0])
#         child = parent1.copy()
#         child[point:] = parent2[point:]
#         return child
    
#     def arithmetic_crossover(self, parent1, parent2):
#         """Arithmetic crossover"""
#         alpha = np.random.rand()
#         return alpha * parent1 + (1 - alpha) * parent2
    
    # def crossover(self, parent1, parent2):
    #     """Perform crossover based on method"""
    #     if self.crossover_method == 'uniform':
    #         return self.uniform_crossover(parent1, parent2)
    #     elif self.crossover_method == 'single_point':
    #         return self.single_point_crossover(parent1, parent2)
    #     elif self.crossover_method == 'arithmetic':
    #         return self.arithmetic_crossover(parent1, parent2)
    #     else:
    #         raise ValueError(f"Unknown crossover method: {self.crossover_method}")
    
#     def gaussian_mutation(self, individual):
#         """Gaussian mutation"""
#         mutation_mask = np.random.rand(*individual.shape) < self.mutation_rate
#         mutations = np.random.randn(*individual.shape) * 0.1
#         individual = individual.copy()
#         individual[mutation_mask] += mutations[mutation_mask]
#         return individual
    
#     def uniform_mutation(self, individual):
#         """Uniform mutation"""
#         mutation_mask = np.random.rand(*individual.shape) < self.mutation_rate
#         mutations = np.random.uniform(-0.2, 0.2, individual.shape)
#         individual = individual.copy()
#         individual[mutation_mask] += mutations[mutation_mask]
#         return individual
    
#     def mutate(self, individual):
#         """Mutate based on method"""
#         if self.mutation_method == 'gaussian':
#             return self.gaussian_mutation(individual)
#         elif self.mutation_method == 'uniform':
#             return self.uniform_mutation(individual)
#         else:
#             raise ValueError(f"Unknown mutation method: {self.mutation_method}")
    
#     def calculate_diversity(self, population):
#         """Calculate population diversity"""
#         if len(population) < 2:
#             return 0.0
        
#         total_distance = 0
#         count = 0
#         for i in range(len(population)):
#             for j in range(i + 1, len(population)):
#                 distance = np.linalg.norm(population[i] - population[j])
#                 total_distance += distance
#                 count += 1
        
#         return total_distance / count if count > 0 else 0.0
    
#     def check_convergence(self, fitness_history, window=10, threshold=1e-4):
#         """Check if algorithm has converged"""
#         if len(fitness_history) < window:
#             return False
        
#         recent_fitness = fitness_history[-window:]
#         variance = np.var(recent_fitness)
#         return variance < threshold
    
#     def run(self):
#         """Run the evolutionary algorithm"""
#         start_time = time.time()
#         self.logger.info(f"Starting EA with config: pop_size={self.population_size}, "
#                         f"generations={self.generations}, selection={self.selection_method}, "
#                         f"crossover={self.crossover_method}, mutation={self.mutation_method}")
        
#         # Initialize population
#         population = [self.generate_individual() for _ in range(self.population_size)]
        
#         for gen in range(self.generations):
#             # Evaluate fitness for all individuals
#             fitness_scores = [self.evaluate(ind) for ind in population]
            
#             # Calculate statistics
#             best_fitness = min(fitness_scores)
#             avg_fitness = np.mean(fitness_scores)
#             diversity = self.calculate_diversity(population)
            
#             # Store statistics
#             self.stats['fitness_history'].append(best_fitness)
#             self.stats['avg_fitness_history'].append(avg_fitness)
#             self.stats['diversity_history'].append(diversity)
            
#             # Check convergence
#             if self.check_convergence(self.stats['fitness_history']) and self.stats['convergence_gen'] is None:
#                 self.stats['convergence_gen'] = gen
#                 self.logger.info(f"Convergence detected at generation {gen}")
            
#             # Log progress
#             if gen % 10 == 0 or gen == self.generations - 1:
#                 self.logger.info(f"Gen {gen + 1}/{self.generations}: Best={best_fitness:.4f}, "
#                                f"Avg={avg_fitness:.4f}, Diversity={diversity:.4f}")
            
#             # Sort by fitness
#             sorted_indices = np.argsort(fitness_scores)
#             sorted_population = [population[i] for i in sorted_indices]
#             sorted_fitness = [fitness_scores[i] for i in sorted_indices]
            
#             # Create new population
#             new_population = []
            
#             # Elitism - keep best individual
#             new_population.append(sorted_population[0].copy())
            
#             # Generate rest of population
#             while len(new_population) < self.population_size:
#                 parent1 = self.select_parent(population, fitness_scores)
#                 parent2 = self.select_parent(population, fitness_scores)
#                 child = self.crossover(parent1, parent2)
#                 child = self.mutate(child)
#                 new_population.append(child)
            
#             population = new_population
        
#         # Final evaluation
#         final_fitness = [self.evaluate(ind) for ind in population]
#         best_idx = np.argmin(final_fitness)
        
#         self.stats['execution_time'] = time.time() - start_time
#         self.logger.info(f"Evolution complete! Final best RMSD: {final_fitness[best_idx]:.4f}")
#         self.logger.info(f"Execution time: {self.stats['execution_time']:.2f} seconds")
        
#         return population[best_idx], self.stats
    
#     def plot_results(self, save_path=None):
#         """Plot evolution results"""
#         fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
        
#         # Fitness evolution
#         ax1.plot(self.stats['fitness_history'], label='Best Fitness', color='blue')
#         ax1.plot(self.stats['avg_fitness_history'], label='Average Fitness', color='red', alpha=0.7)
#         ax1.set_xlabel('Generation')
#         ax1.set_ylabel('RMSD')
#         ax1.set_title('Fitness Evolution')
#         ax1.legend()
#         ax1.grid(True)
        
#         # Diversity
#         ax2.plot(self.stats['diversity_history'], color='green')
#         ax2.set_xlabel('Generation')
#         ax2.set_ylabel('Population Diversity')
#         ax2.set_title('Population Diversity')
#         ax2.grid(True)
        
#         # Convergence analysis
#         improvement = np.diff(self.stats['fitness_history'])
#         ax3.plot(improvement, color='purple')
#         ax3.set_xlabel('Generation')
#         ax3.set_ylabel('Fitness Improvement')
#         ax3.set_title('Fitness Improvement per Generation')
#         ax3.grid(True)
        
#         # Selection pressure
#         if len(self.stats['fitness_history']) > 1:
#             selection_pressure = []
#             for i in range(1, len(self.stats['fitness_history'])):
#                 if self.stats['avg_fitness_history'][i] != 0:
#                     pressure = self.stats['fitness_history'][i] / self.stats['avg_fitness_history'][i]
#                     selection_pressure.append(pressure)
#                 else:
#                     selection_pressure.append(1.0)
#             ax4.plot(selection_pressure, color='orange')
#             ax4.set_xlabel('Generation')
#             ax4.set_ylabel('Selection Pressure (Best/Avg)')
#             ax4.set_title('Selection Pressure')
#             ax4.grid(True)
        
#         plt.tight_layout()
#         if save_path:
#             plt.savefig(save_path, dpi=300, bbox_inches='tight')
#         plt.show()
    
#     def get_summary_stats(self):
#         """Get summary statistics"""
#         return {
#             'final_best_fitness': self.stats['fitness_history'][-1] if self.stats['fitness_history'] else None,
#             'initial_fitness': self.stats['fitness_history'][0] if self.stats['fitness_history'] else None,
#             'improvement': (self.stats['fitness_history'][0] - self.stats['fitness_history'][-1]) if len(self.stats['fitness_history']) > 1 else 0,
#             'convergence_generation': self.stats['convergence_gen'],
#             'execution_time': self.stats['execution_time'],
#             'final_diversity': self.stats['diversity_history'][-1] if self.stats['diversity_history'] else None,
#             'avg_diversity': np.mean(self.stats['diversity_history']) if self.stats['diversity_history'] else None
#         }
    
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
