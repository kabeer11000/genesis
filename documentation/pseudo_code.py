import numpy as np
import torch
import torch.nn as nn
import random
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from enum import Enum
import copy

# ============================================================================
# DATA STRUCTURES AND ENUMS
# ============================================================================

class FunctionType(Enum):
    LINEAR = "linear"
    NONLINEAR = "nonlinear"
    ATTENTION = "attention"
    CONV = "conv"
    LSTM = "lstm"
    TRANSFORMER = "transformer"

@dataclass
class Function:
    type: FunctionType
    weights: torch.Tensor
    input_refs: List[int]  # indices of functions this uses as input
    coefficient: float     # 0-1 weight in final ensemble
    network: Optional[nn.Module] = None

@dataclass
class Individual:
    predefined_functions: List[Function]
    evolved_functions: List[Function]
    fitness: float = 0.0

@dataclass
class ProteinData:
    sequence: str
    structure: torch.Tensor  # 3D coordinates
    features: torch.Tensor

# ============================================================================
# PREDEFINED FUNCTION ARCHITECTURES
# ============================================================================

class PhiPsiPredictor(nn.Module):
    """Predicts backbone dihedral angles (phi, psi)"""
    def __init__(self, input_dim: int = 20):
        super().__init__()
        self.lstm = nn.LSTM(input_dim, 128, batch_first=True, bidirectional=True)
        self.phi_head = nn.Linear(256, 1)
        self.psi_head = nn.Linear(256, 1)
    
    def forward(self, sequence_features):
        lstm_out, _ = self.lstm(sequence_features)
        phi_angles = self.phi_head(lstm_out)
        psi_angles = self.psi_head(lstm_out)
        return torch.cat([phi_angles, psi_angles], dim=-1)

class DistanceMatrixPredictor(nn.Module):
    """Predicts inter-residue distance matrix"""
    def __init__(self, input_dim: int = 20):
        super().__init__()
        self.feature_net = nn.Sequential(
            nn.Linear(input_dim, 64),
            nn.ReLU(),
            nn.Linear(64, 32)
        )
        self.distance_head = nn.Linear(64, 1)  # 64 = 32*2 for pairwise
    
    def forward(self, sequence_features):
        batch_size, seq_len, _ = sequence_features.shape
        features = self.feature_net(sequence_features)
        
        # Create pairwise features
        left = features.unsqueeze(2).expand(-1, -1, seq_len, -1)
        right = features.unsqueeze(1).expand(-1, seq_len, -1, -1)
        pairwise = torch.cat([left, right], dim=-1)
        
        distances = self.distance_head(pairwise).squeeze(-1)
        return distances

class SecondaryStructurePredictor(nn.Module):
    """Predicts secondary structure probabilities"""
    def __init__(self, input_dim: int = 20):
        super().__init__()
        self.transformer = nn.TransformerEncoder(
            nn.TransformerEncoderLayer(input_dim, nhead=4, dim_feedforward=128),
            num_layers=2
        )
        self.ss_head = nn.Linear(input_dim, 3)  # helix, sheet, loop
    
    def forward(self, sequence_features):
        transformer_out = self.transformer(sequence_features.transpose(0, 1)).transpose(0, 1)
        ss_probs = torch.softmax(self.ss_head(transformer_out), dim=-1)
        return ss_probs

# ============================================================================
# DYNAMIC FUNCTION GENERATORS
# ============================================================================

class DynamicLinear(nn.Module):
    def __init__(self, input_dims: List[int], output_dim: int = 64):
        super().__init__()
        total_input = sum(input_dims) if input_dims else 64
        self.linear = nn.Linear(total_input, output_dim)
    
    def forward(self, inputs: List[torch.Tensor]):
        if inputs:
            concatenated = torch.cat([inp.flatten(start_dim=1) for inp in inputs], dim=-1)
        else:
            concatenated = torch.zeros(1, 64)
        return self.linear(concatenated)

class DynamicAttention(nn.Module):
    def __init__(self, input_dims: List[int], output_dim: int = 64):
        super().__init__()
        self.attention = nn.MultiheadAttention(output_dim, num_heads=4)
        self.projection = nn.Linear(sum(input_dims) if input_dims else 64, output_dim)
    
    def forward(self, inputs: List[torch.Tensor]):
        if inputs:
            concatenated = torch.cat([inp.flatten(start_dim=1) for inp in inputs], dim=-1)
        else:
            concatenated = torch.zeros(1, 64)
        
        projected = self.projection(concatenated).unsqueeze(0)
        attn_out, _ = self.attention(projected, projected, projected)
        return attn_out.squeeze(0)

# ============================================================================
# NEUROEVOLUTIONARY SYSTEM
# ============================================================================

class ProteinStructurePredictor:
    def __init__(self, population_size: int = 50, num_predefined: int = 3, 
                 num_evolved: int = 7, max_generations: int = 100):
        self.population_size = population_size
        self.k = num_predefined  # number of predefined functions
        self.n = num_predefined + num_evolved  # total functions
        self.max_generations = max_generations
        self.elite_size = 5
        self.pso_iterations = 20
        
    def generate_random_function(self, function_index: int) -> Function:
        """Generate a random evolved function"""
        func_type = random.choice(list(FunctionType))
        
        # Can reference any previous function
        available_refs = list(range(function_index))
        input_refs = random.sample(available_refs, 
                                 min(random.randint(0, 3), len(available_refs)))
        
        # Create network based on type
        if func_type == FunctionType.LINEAR:
            network = DynamicLinear([64] * len(input_refs))
        elif func_type == FunctionType.ATTENTION:
            network = DynamicAttention([64] * len(input_refs))
        else:
            network = DynamicLinear([64] * len(input_refs))  # fallback
        
        return Function(
            type=func_type,
            weights=torch.randn(100),  # generic weight tensor
            input_refs=input_refs,
            coefficient=random.random(),
            network=network
        )
    
    def initialize_population(self) -> List[Individual]:
        """Initialize population with predefined + evolved functions"""
        population = []
        
        for _ in range(self.population_size):
            # Initialize predefined functions
            predefined = [
                Function(FunctionType.LINEAR, torch.randn(100), [], random.random(), 
                        PhiPsiPredictor()),
                Function(FunctionType.LINEAR, torch.randn(100), [], random.random(), 
                        DistanceMatrixPredictor()),
                Function(FunctionType.LINEAR, torch.randn(100), [], random.random(), 
                        SecondaryStructurePredictor())
            ]
            
            # Generate evolved functions
            evolved = []
            for i in range(self.k, self.n):
                evolved.append(self.generate_random_function(i))
            
            population.append(Individual(predefined, evolved))
        
        return population
    
    def evaluate_fitness(self, individual: Individual, training_data: List[ProteinData]) -> float:
        """Evaluate fitness of an individual on training data"""
        total_error = 0.0
        
        for protein in training_data[:10]:  # sample for speed
            try:
                predicted = self.predict_structure(individual, protein.sequence, protein.features)
                true_structure = protein.structure
                
                # Multi-objective fitness
                rmsd_error = self.calculate_rmsd(predicted, true_structure)
                contact_accuracy = self.calculate_contact_accuracy(predicted, true_structure)
                
                error = 0.7 * rmsd_error + 0.3 * (1 - contact_accuracy)
                total_error += error
                
            except Exception as e:
                total_error += 100  # penalty for failed predictions
        
        return 1.0 / (1.0 + total_error)
    
    def predict_structure(self, individual: Individual, sequence: str, 
                         sequence_features: torch.Tensor) -> torch.Tensor:
        """Predict 3D structure using individual's function composition"""
        function_outputs = {}
        
        # Execute predefined functions
        for i, func in enumerate(individual.predefined_functions):
            if i == 0:  # PhiPsi predictor
                function_outputs[i] = func.network(sequence_features.unsqueeze(0))
            elif i == 1:  # Distance matrix predictor
                function_outputs[i] = func.network(sequence_features.unsqueeze(0))
            elif i == 2:  # Secondary structure predictor
                function_outputs[i] = func.network(sequence_features.unsqueeze(0))
        
        # Execute evolved functions
        for i, func in enumerate(individual.evolved_functions, start=self.k):
            # Gather inputs from referenced functions
            inputs = []
            for ref_idx in func.input_refs:
                if ref_idx in function_outputs:
                    inputs.append(function_outputs[ref_idx])
            
            # Execute function
            if inputs:
                function_outputs[i] = func.network(inputs)
            else:
                # Fallback if no valid inputs
                function_outputs[i] = torch.zeros(1, 64)
        
        # Combine outputs with coefficients
        combined_output = self.weighted_combination(function_outputs, individual)
        
        # Convert to 3D coordinates (simplified)
        seq_len = len(sequence)
        structure_3d = combined_output[:, :seq_len*3].reshape(-1, seq_len, 3)
        
        return structure_3d
    
    def weighted_combination(self, function_outputs: Dict[int, torch.Tensor], 
                           individual: Individual) -> torch.Tensor:
        """Combine function outputs with learned coefficients"""
        combined = torch.zeros(1, 1000)  # large enough output
        
        # Combine predefined functions
        for i, func in enumerate(individual.predefined_functions):
            if i in function_outputs:
                output = function_outputs[i].flatten()
                size = min(output.size(0), combined.size(1))
                combined[0, :size] += func.coefficient * output[:size]
        
        # Combine evolved functions
        for i, func in enumerate(individual.evolved_functions, start=self.k):
            if i in function_outputs:
                output = function_outputs[i].flatten()
                size = min(output.size(0), combined.size(1))
                combined[0, :size] += func.coefficient * output[:size]
        
        return combined
    
    def crossover(self, parent1: Individual, parent2: Individual) -> Individual:
        """Create offspring through crossover"""
        child = Individual([], [])
        
        # Crossover predefined functions (blend weights)
        for i in range(self.k):
            p1_func = parent1.predefined_functions[i]
            p2_func = parent2.predefined_functions[i]
            
            # Blend coefficients
            alpha = random.random()
            new_coeff = alpha * p1_func.coefficient + (1-alpha) * p2_func.coefficient
            
            child_func = copy.deepcopy(p1_func)
            child_func.coefficient = new_coeff
            child.predefined_functions.append(child_func)
        
        # Crossover evolved functions (structure + weights)
        for i in range(len(parent1.evolved_functions)):
            if random.random() < 0.5:
                child.evolved_functions.append(copy.deepcopy(parent1.evolved_functions[i]))
            else:
                child.evolved_functions.append(copy.deepcopy(parent2.evolved_functions[i]))
        
        return child
    
    def mutate(self, individual: Individual) -> Individual:
        """Mutate an individual"""
        # Mutate predefined function coefficients
        for func in individual.predefined_functions:
            if random.random() < 0.1:
                func.coefficient += random.gauss(0, 0.1)
                func.coefficient = max(0, min(1, func.coefficient))
        
        # Mutate evolved functions
        for i, func in enumerate(individual.evolved_functions):
            if random.random() < 0.2:
                # Structure mutation: change input references
                if random.random() < 0.5:
                    available_refs = list(range(self.k + i))
                    func.input_refs = random.sample(available_refs, 
                                                  min(random.randint(0, 3), len(available_refs)))
                
                # Coefficient mutation
                func.coefficient += random.gauss(0, 0.1)
                func.coefficient = max(0, min(1, func.coefficient))
        
        return individual
    
    def tournament_selection(self, population: List[Individual], k: int = 3) -> Individual:
        """Tournament selection"""
        tournament = random.sample(population, min(k, len(population)))
        return max(tournament, key=lambda x: x.fitness)
    
    def train(self, training_data: List[ProteinData]) -> Individual:
        """Main training loop"""
        population = self.initialize_population()
        
        for generation in range(self.max_generations):
            # Evaluate fitness
            for individual in population:
                individual.fitness = self.evaluate_fitness(individual, training_data)
            
            # Sort by fitness
            population.sort(key=lambda x: x.fitness, reverse=True)
            
            print(f"Generation {generation}: Best fitness = {population[0].fitness:.4f}")
            
            # Create new population
            new_population = []
            
            # Elite selection
            new_population.extend(population[:self.elite_size])
            
            # Generate offspring
            while len(new_population) < self.population_size:
                parent1 = self.tournament_selection(population)
                parent2 = self.tournament_selection(population)
                child = self.crossover(parent1, parent2)
                child = self.mutate(child)
                new_population.append(child)
            
            population = new_population
        
        # Return best individual
        for individual in population:
            individual.fitness = self.evaluate_fitness(individual, training_data)
        
        return max(population, key=lambda x: x.fitness)
    
    def predict_novel_protein(self, trained_individual: Individual, 
                            novel_sequence: str, sequence_features: torch.Tensor) -> Dict:
        """Predict structure of novel protein"""
        predicted_structure = self.predict_structure(trained_individual, novel_sequence, sequence_features)
        
        # Estimate confidence (simplified)
        confidence = min(trained_individual.fitness, 1.0)
        
        return {
            'structure': predicted_structure,
            'confidence': confidence,
            'sequence': novel_sequence
        }
    
    # Helper methods for fitness calculation
    def calculate_rmsd(self, pred: torch.Tensor, true: torch.Tensor) -> float:
        """Calculate RMSD between predicted and true structures"""
        if pred.shape != true.shape:
            return 100.0  # high penalty for shape mismatch
        
        diff = pred - true
        rmsd = torch.sqrt(torch.mean(diff ** 2))
        return rmsd.item()
    
    def calculate_contact_accuracy(self, pred: torch.Tensor, true: torch.Tensor) -> float:
        """Calculate contact map accuracy (simplified)"""
        # This is a simplified version - in practice you'd compute actual contact maps
        return random.random()  # placeholder

# ============================================================================
# USAGE EXAMPLE
# ============================================================================

def create_dummy_data(num_proteins: int = 100) -> List[ProteinData]:
    """Create dummy training data for testing"""
    data = []
    for i in range(num_proteins):
        seq_len = random.randint(50, 200)
        sequence = ''.join(random.choices('ACDEFGHIKLMNPQRSTVWY', k=seq_len))
        features = torch.randn(seq_len, 20)  # 20D amino acid features
        structure = torch.randn(seq_len, 3)  # 3D coordinates
        
        data.append(ProteinData(sequence, structure, features))
    
    return data

def main():
    print("Initializing Neuroevolutionary Protein Structure Predictor...")
    predictor = ProteinStructurePredictor(
        population_size=20,  # smaller for demo
        num_predefined=3,
        num_evolved=5,
        max_generations=10  # fewer generations for demo
    )
    
    print("Creating training data...")
    training_data = create_dummy_data(50)
    
    print("Training the system...")
    best_individual = predictor.train(training_data)
    
    print(f"Training completed! Best fitness: {best_individual.fitness:.4f}")
    
    # Predict novel protein
    novel_sequence = "MKLLVLLLLVLVAAAATAAQGDDFDKFLTKYAELKSIEESEDLKL"
    novel_features = torch.randn(len(novel_sequence), 20)
    
    prediction = predictor.predict_novel_protein(best_individual, novel_sequence, novel_features)
    
    print(f"Novel protein prediction completed!")
    print(f"Sequence length: {len(novel_sequence)}")
    print(f"Predicted structure shape: {prediction['structure'].shape}")
    print(f"Confidence: {prediction['confidence']:.4f}")

if __name__ == "__main__":
    main()