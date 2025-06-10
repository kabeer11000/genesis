from ea import EvolutionaryAlgorithm
# from pso import ParticleSwarmOptimizer
from utils import load_proteins

# Load protein dataset
protein_paths = ['data/1A8O.pdb', 'data/1CRN.pdb', 'data/1FME.pdb']  # Add your paths here
proteins = load_proteins(protein_paths)

# Run EA
ea = EvolutionaryAlgorithm(proteins)
ea.run()

# Run PSO
# pso = ParticleSwarmOptimizer(proteins)
# pso.run()
