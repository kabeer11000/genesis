#!/bin/bash

# Create Genesis Project Setup Script
# This script creates all the necessary files for the Genesis project

echo "Creating Genesis Project structure..."

# Create data directory
mkdir -p data
echo "✓ Created data/ directory"

# Create main.py
cat > main.py << 'EOF'
from ea import EvolutionaryAlgorithm
from utils import load_proteins

# Load protein dataset
protein_paths = ['data/1A8O.pdb', 'data/1CRN.pdb', 'data/1FME.pdb']  # Add your paths here
proteins = load_proteins(protein_paths)

# Run EA
ea = EvolutionaryAlgorithm(proteins)
ea.run()
EOF
echo "✓ Created main.py"

# Create utils.py
cat > utils.py << 'EOF'
from Bio.PDB import PDBParser, PPBuilder
import numpy as np


def load_proteins(paths):
    parser = PDBParser(QUIET=True)
    proteins = []
    for path in paths:
        structure = parser.get_structure('protein', path)
        model = structure[0]
        chain = next(model.get_chains())
        
        # Get CA atoms and coordinates
        ca_atoms = []
        coords = []
        for residue in chain:
            if 'CA' in residue:
                ca_atoms.append(residue['CA'])
                coords.append(residue['CA'].get_coord())
        
        # Get sequence using PPBuilder for single-letter codes
        ppb = PPBuilder()
        peptides = ppb.build_peptides(chain)
        seq = ""
        for peptide in peptides:
            seq += str(peptide.get_sequence())
        
        proteins.append({
            'structure': structure,
            'ca_atoms': ca_atoms,
            'coords': np.array(coords),
            'sequence': seq
        })
    return proteins


def rmsd(pred_coords, true_coords):
    """Calculate RMSD between two sets of coordinates"""
    if pred_coords.shape != true_coords.shape:
        raise ValueError("Coordinate arrays must have the same shape")
    
    # Center both structures
    pred_centered = pred_coords - np.mean(pred_coords, axis=0)
    true_centered = true_coords - np.mean(true_coords, axis=0)
    
    # Calculate RMSD
    diff = pred_centered - true_centered
    rmsd_value = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    return rmsd_value
EOF
echo "✓ Created utils.py"

# Create features.py
cat > features.py << 'EOF'
import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def feature_hydrophobicity(protein):
    """Extract hydrophobicity-based features"""
    seq = protein['sequence']
    if not seq:
        return np.zeros(10)
    
    try:
        analysis = ProteinAnalysis(seq)
        # Create a feature vector with hydrophobicity-related metrics
        features = [
            analysis.gravy(),  # Grand average of hydropathy
            analysis.aromaticity(),
            analysis.instability_index(),
            analysis.isoelectric_point(),
            analysis.charge_at_pH(7.0),
        ]
        # Pad to 10 features
        features.extend([0.0] * (10 - len(features)))
        return np.array(features[:10])
    except:
        return np.zeros(10)


def feature_molecular_weight(protein):
    """Extract molecular weight and composition features"""
    seq = protein['sequence']
    if not seq:
        return np.zeros(10)
    
    try:
        analysis = ProteinAnalysis(seq)
        # Get amino acid percentages for common residues
        aa_percent = analysis.get_amino_acids_percent()
        features = [
            analysis.molecular_weight() / 10000.0,  # Normalize
            aa_percent.get('A', 0.0),  # Alanine
            aa_percent.get('L', 0.0),  # Leucine
            aa_percent.get('G', 0.0),  # Glycine
            aa_percent.get('V', 0.0),  # Valine
            aa_percent.get('E', 0.0),  # Glutamate
            aa_percent.get('K', 0.0),  # Lysine
            aa_percent.get('D', 0.0),  # Aspartate
            aa_percent.get('S', 0.0),  # Serine
            aa_percent.get('T', 0.0),  # Threonine
        ]
        return np.array(features[:10])
    except:
        return np.zeros(10)


def feature_secondary_structure(protein):
    """Estimate secondary structure propensities"""
    seq = protein['sequence']
    if not seq:
        return np.zeros(10)
    
    try:
        analysis = ProteinAnalysis(seq)
        # Get secondary structure fractions
        helix, turn, sheet = analysis.secondary_structure_fraction()
        
        # Additional structural features
        features = [
            helix,
            turn,
            sheet,
            len(seq) / 100.0,  # Normalized length
            analysis.molar_extinction_coefficient()[0] / 10000.0,
            analysis.molar_extinction_coefficient()[1] / 10000.0,
            0.0, 0.0, 0.0, 0.0  # Padding
        ]
        return np.array(features[:10])
    except:
        return np.zeros(10)


FEATURE_FUNCTIONS = [
    feature_hydrophobicity,
    feature_molecular_weight,
    feature_secondary_structure
]
EOF
echo "✓ Created features.py"

# Create ea.py
cat > ea.py << 'EOF'
import numpy as np
from features import FEATURE_FUNCTIONS
from utils import rmsd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


class EvolutionaryAlgorithm:
    def __init__(self, proteins, population_size=30, generations=50, mutation_rate=0.1):
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
EOF
echo "✓ Created ea.py"

# Create requirements.txt
cat > requirements.txt << 'EOF'
biopython>=1.79
numpy>=1.21.0
scikit-learn>=1.0.0
tqdm>=4.62.0
EOF
echo "✓ Created requirements.txt"

# Create README.md
cat > README.md << 'EOF'
# Genesis Project - Protein Structure Prediction using Evolutionary Algorithms

## Setup

1. Create a virtual environment:
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

3. Add your PDB files to the `data/` directory

4. Update the file paths in `main.py` to point to your PDB files

5. Run the program:
   ```bash
   python main.py
   ```

## Note

This is a demonstration of evolutionary algorithms applied to feature-to-coordinate mapping. 
For actual protein structure prediction, use specialized tools like AlphaFold or RoseTTAFold.
EOF
echo "✓ Created README.md"

echo ""
echo "Genesis Project structure created successfully!"
echo ""
echo "Directory structure:"
echo "."
echo "├── main.py"
echo "├── ea.py"
echo "├── features.py"
echo "├── utils.py"
echo "├── requirements.txt"
echo "├── README.md"
echo "└── data/"
echo "    └── [place your .pdb files here]"
echo ""
echo "Next steps:"
echo "1. Create a virtual environment: python -m venv .venv"
echo "2. Activate it: source .venv/bin/activate"
echo "3. Install dependencies: pip install -r requirements.txt"
echo "4. Add PDB files to the data/ directory"
echo "5. Run: python main.py"