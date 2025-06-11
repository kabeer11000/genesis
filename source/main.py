# from ea import EvolutionaryAlgorithm
# from utils import load_proteins
# import matplotlib.pyplot as plt

# def run_basic_example():
#     """Run a basic example with the enhanced EA"""
#     print("Loading proteins...")
    
#     # Try to load real protein files, otherwise create dummy data
#     try:
#         protein_paths = ['data/1A8O.pdb', 'data/1CRN.pdb', 'data/1FME.pdb']
#         proteins = load_proteins(protein_paths)
#         print(f"Loaded {len(proteins)} real proteins")
#     except:
#         print("Real protein files not found, creating dummy data...")
#         import numpy as np
        
#         # Create dummy protein data for testing
#         dummy_proteins = []
#         for i in range(3):
#             n_residues = np.random.randint(50, 100)
#             coords = np.random.randn(n_residues, 3) * 10
#             sequence = ''.join(np.random.choice(list('ACDEFGHIKLMNPQRSTVWY'), size=n_residues))
            
#             dummy_proteins.append({
#                 'coords': coords,
#                 'sequence': sequence,
#                 'ca_atoms': None,
#                 'structure': None
#             })
        
#         proteins = dummy_proteins
#         print(f"Created {len(proteins)} dummy proteins for testing")
    
#     # Run EA with different configurations
#     configurations = [
#         {
#             'name': 'Tournament Selection',
#             'selection_method': 'tournament',
#             'crossover_method': 'uniform',
#             'mutation_method': 'gaussian',
#             'population_size': 50,
#             'generations': 25
#         },
#         {
#             'name': 'Binary Tournament',
#             'selection_method': 'binary_tournament',
#             'crossover_method': 'arithmetic',
#             'mutation_method': 'gaussian',
#             'population_size': 50,
#             'generations': 25
#         }
#     ]
    
#     results = []
    
#     for config in configurations:
#         print(f"\n{'='*50}")
#         print(f"Running: {config['name']}")
#         print(f"{'='*50}")
        
#         # Create and run EA
#         ea = EvolutionaryAlgorithm(
#             proteins,
#             population_size=config['population_size'],
#             generations=config['generations'],
#             selection_method=config['selection_method'],
#             crossover_method=config['crossover_method'],
#             mutation_method=config['mutation_method'],
#             log_level='INFO'
#         )
        
#         best_individual, stats = ea.run()
#         summary = ea.get_summary_stats()
        
#         print(f"\nResults for {config['name']}:")
#         print(f"  Final Best RMSD: {summary['final_best_fitness']:.4f}")
#         print(f"  Improvement: {summary['improvement']:.4f}")
#         print(f"  Execution Time: {summary['execution_time']:.2f}s")
#         if summary['convergence_generation']:
#             print(f"  Converged at generation: {summary['convergence_generation']}")
        
#         # Plot results
#         ea.plot_results(save_path=f"results_{config['name'].replace(' ', '_').lower()}.png")
        
#         results.append({
#             'config': config['name'],
#             'stats': stats,
#             'summary': summary
#         })
    
#     # Compare results
#     print(f"\n{'='*50}")
#     print("COMPARISON SUMMARY")
#     print(f"{'='*50}")
    
#     for result in results:
#         print(f"{result['config']}: Final RMSD = {result['summary']['final_best_fitness']:.4f}")
    
#     # Plot comparison
#     plt.figure(figsize=(12, 6))
#     for result in results:
#         plt.plot(result['stats']['fitness_history'], 
#                 label=result['config'], linewidth=2)
    
#     plt.xlabel('Generation')
#     plt.ylabel('Best RMSD')
#     plt.title('EA Configuration Comparison')
#     plt.legend()
#     plt.grid(True, alpha=0.3)
#     plt.tight_layout()
#     plt.savefig('ea_comparison.png', dpi=300, bbox_inches='tight')
#     plt.show()
    
#     print("\nExample completed! Check the generated plots and log files.")

# if __name__ == "__main__":
#     run_basic_example()

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
