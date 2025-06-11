import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
try:
    import seaborn as sns
except ImportError:
    print("Seaborn not available, using matplotlib only")
    sns = None
from ea import EvolutionaryAlgorithm
from utils import load_proteins
import time
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Set style for better plots
try:
    plt.style.use('seaborn-v0_8')
except:
    try:
        plt.style.use('seaborn')
    except:
        plt.style.use('default')

if sns is not None:
    sns.set_palette("husl")

def create_test_configs():
    """Create different EA configurations to test"""
    configs = [
        # Basic configurations
        {'name': 'Tournament_Uniform_Gaussian', 'selection_method': 'tournament', 'crossover_method': 'uniform', 'mutation_method': 'gaussian'},
        {'name': 'Binary_Tournament_Uniform_Gaussian', 'selection_method': 'binary_tournament', 'crossover_method': 'uniform', 'mutation_method': 'gaussian'},
        {'name': 'Roulette_Uniform_Gaussian', 'selection_method': 'roulette', 'crossover_method': 'uniform', 'mutation_method': 'gaussian'},
        {'name': 'Rank_Uniform_Gaussian', 'selection_method': 'rank', 'crossover_method': 'uniform', 'mutation_method': 'gaussian'},
        
        # Crossover method comparisons
        {'name': 'Tournament_SinglePoint_Gaussian', 'selection_method': 'tournament', 'crossover_method': 'single_point', 'mutation_method': 'gaussian'},
        {'name': 'Tournament_Arithmetic_Gaussian', 'selection_method': 'tournament', 'crossover_method': 'arithmetic', 'mutation_method': 'gaussian'},
        
        # Mutation method comparisons
        {'name': 'Tournament_Uniform_Uniform_Mut', 'selection_method': 'tournament', 'crossover_method': 'uniform', 'mutation_method': 'uniform'},
        
        # Population size variations
        {'name': 'Small_Pop_50', 'selection_method': 'tournament', 'crossover_method': 'uniform', 'mutation_method': 'gaussian', 'population_size': 50},
        {'name': 'Large_Pop_200', 'selection_method': 'tournament', 'crossover_method': 'uniform', 'mutation_method': 'gaussian', 'population_size': 200},
        
        # Mutation rate variations
        {'name': 'Low_Mutation_0.05', 'selection_method': 'tournament', 'crossover_method': 'uniform', 'mutation_method': 'gaussian', 'mutation_rate': 0.05},
        {'name': 'High_Mutation_0.2', 'selection_method': 'tournament', 'crossover_method': 'uniform', 'mutation_method': 'gaussian', 'mutation_rate': 0.2},
    ]
    
    return configs

def run_single_test(config, proteins, generations=30, trials=3):
    """Run a single EA configuration multiple times"""
    print(f"\n{'='*50}")
    print(f"Testing: {config['name']}")
    print(f"{'='*50}")
    
    results = []
    all_fitness_histories = []
    
    for trial in range(trials):
        print(f"\nTrial {trial + 1}/{trials}")
        
        # Set default parameters
        params = {
            'population_size': 100,
            'generations': generations,
            'mutation_rate': 0.1,
            'log_level': 'WARNING'  # Reduce logging during testing
        }
        params.update(config)
        params.pop('name', None)  # Remove name from parameters
        
        # Run EA
        ea = EvolutionaryAlgorithm(proteins, **params)
        best_individual, stats = ea.run()
        
        # Store results
        summary = ea.get_summary_stats()
        summary['trial'] = trial + 1
        summary['config'] = config['name']
        results.append(summary)
        all_fitness_histories.append(stats['fitness_history'])
    
    return results, all_fitness_histories

def analyze_results(all_results):
    """Analyze and summarize results from all tests"""
    df = pd.DataFrame(all_results)
    
    # Group by configuration
    summary_stats = df.groupby('config').agg({
        'final_best_fitness': ['mean', 'std', 'min', 'max'],
        'improvement': ['mean', 'std'],
        'convergence_generation': ['mean', 'std'],
        'execution_time': ['mean', 'std'],
        'final_diversity': ['mean', 'std'],
        'avg_diversity': ['mean', 'std']
    }).round(4)
    
    return df, summary_stats

def plot_comparison_results(all_results, all_histories, save_dir='results'):
    """Create comprehensive comparison plots"""
    os.makedirs(save_dir, exist_ok=True)
    
    df = pd.DataFrame(all_results)
    
    # 1. Final fitness comparison
    plt.figure(figsize=(14, 8))
    if sns is not None:
        try:
            sns.boxplot(data=df, x='config', y='final_best_fitness')
        except:
            # Fallback to matplotlib if seaborn fails
            configs = df['config'].unique()
            data_by_config = [df[df['config'] == config]['final_best_fitness'].values for config in configs]
            plt.boxplot(data_by_config, labels=configs)
    else:
        configs = df['config'].unique()
        data_by_config = [df[df['config'] == config]['final_best_fitness'].values for config in configs]
        plt.boxplot(data_by_config, labels=configs)
    
    plt.xticks(rotation=45, ha='right')
    plt.title('Final Best Fitness Comparison Across Configurations')
    plt.ylabel('Final RMSD')
    plt.tight_layout()
    plt.savefig(f'{save_dir}/final_fitness_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. Improvement comparison
    plt.figure(figsize=(14, 8))
    if sns is not None:
        try:
            sns.boxplot(data=df, x='config', y='improvement')
        except:
            configs = df['config'].unique()
            data_by_config = [df[df['config'] == config]['improvement'].values for config in configs]
            plt.boxplot(data_by_config, labels=configs)
    else:
        configs = df['config'].unique()
        data_by_config = [df[df['config'] == config]['improvement'].values for config in configs]
        plt.boxplot(data_by_config, labels=configs)
    
    plt.xticks(rotation=45, ha='right')
    plt.title('Fitness Improvement Comparison')
    plt.ylabel('Improvement (Initial - Final RMSD)')
    plt.tight_layout()
    plt.savefig(f'{save_dir}/improvement_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 3. Execution time comparison
    plt.figure(figsize=(14, 8))
    if sns is not None:
        try:
            sns.boxplot(data=df, x='config', y='execution_time')
        except:
            configs = df['config'].unique()
            data_by_config = [df[df['config'] == config]['execution_time'].values for config in configs]
            plt.boxplot(data_by_config, labels=configs)
    else:
        configs = df['config'].unique()
        data_by_config = [df[df['config'] == config]['execution_time'].values for config in configs]
        plt.boxplot(data_by_config, labels=configs)
    
    plt.xticks(rotation=45, ha='right')
    plt.title('Execution Time Comparison')
    plt.ylabel('Time (seconds)')
    plt.tight_layout()
    plt.savefig(f'{save_dir}/execution_time_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 4. Convergence generation comparison
    plt.figure(figsize=(14, 8))
    df_conv = df.dropna(subset=['convergence_generation'])
    if not df_conv.empty:
        if sns is not None:
            try:
                sns.boxplot(data=df_conv, x='config', y='convergence_generation')
            except:
                configs = df_conv['config'].unique()
                data_by_config = [df_conv[df_conv['config'] == config]['convergence_generation'].values for config in configs]
                plt.boxplot(data_by_config, labels=configs)
        else:
            configs = df_conv['config'].unique()
            data_by_config = [df_conv[df_conv['config'] == config]['convergence_generation'].values for config in configs]
            plt.boxplot(data_by_config, labels=configs)
        
        plt.xticks(rotation=45, ha='right')
        plt.title('Convergence Generation Comparison')
        plt.ylabel('Generation')
        plt.tight_layout()
        plt.savefig(f'{save_dir}/convergence_comparison.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    # 5. Fitness evolution curves
    plt.figure(figsize=(16, 10))
    config_names = df['config'].unique()
    colors = plt.cm.tab10(np.linspace(0, 1, len(config_names)))
    
    for i, config in enumerate(config_names):
        config_histories = [hist for j, hist in enumerate(all_histories) 
                          if all_results[j]['config'] == config]
        
        if config_histories:
            # Calculate mean and std across trials
            min_length = min(len(hist) for hist in config_histories)
            trimmed_histories = [hist[:min_length] for hist in config_histories]
            mean_history = np.mean(trimmed_histories, axis=0)
            std_history = np.std(trimmed_histories, axis=0)
            
            generations = range(len(mean_history))
            plt.plot(generations, mean_history, label=config, linewidth=2, color=colors[i])
            plt.fill_between(generations, 
                           mean_history - std_history, 
                           mean_history + std_history, 
                           alpha=0.2, color=colors[i])
    
    plt.xlabel('Generation')
    plt.ylabel('Best RMSD')
    plt.title('Fitness Evolution Comparison (Mean ± Std)')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{save_dir}/fitness_evolution_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 6. Performance heatmap
    plt.figure(figsize=(12, 8))
    metrics = ['final_best_fitness', 'improvement', 'execution_time', 'final_diversity']
    available_metrics = [m for m in metrics if m in df.columns]
    
    if available_metrics:
        heatmap_data = df.groupby('config')[available_metrics].mean()
        
        # Normalize data for better visualization
        heatmap_normalized = (heatmap_data - heatmap_data.min()) / (heatmap_data.max() - heatmap_data.min())
        
        try:
            if sns is not None:
                sns.heatmap(heatmap_normalized.T, annot=True, cmap='RdYlBu_r', 
                            xticklabels=[name.replace('_', '\n') for name in heatmap_data.index],
                            yticklabels=[m.replace('_', ' ').title() for m in available_metrics])
            else:
                # Fallback to matplotlib imshow
                import matplotlib.pyplot as plt
                im = plt.imshow(heatmap_normalized.T, cmap='RdYlBu_r', aspect='auto')
                plt.colorbar(im)
                plt.xticks(range(len(heatmap_data.index)), 
                          [name.replace('_', '\n') for name in heatmap_data.index], rotation=45)
                plt.yticks(range(len(available_metrics)), 
                          [m.replace('_', ' ').title() for m in available_metrics])
        except:
            # Additional fallback to matplotlib imshow
            im = plt.imshow(heatmap_normalized.T, cmap='RdYlBu_r', aspect='auto')
            plt.colorbar(im)
            plt.xticks(range(len(heatmap_data.index)), 
                      [name.replace('_', '\n') for name in heatmap_data.index], rotation=45)
            plt.yticks(range(len(available_metrics)), 
                      [m.replace('_', ' ').title() for m in available_metrics])
        
        plt.title('Performance Metrics Heatmap (Normalized)')
        plt.tight_layout()
        plt.savefig(f'{save_dir}/performance_heatmap.png', dpi=300, bbox_inches='tight')
        plt.show()

def create_summary_tables(df, summary_stats, save_dir='results'):
    """Create and save summary tables"""
    os.makedirs(save_dir, exist_ok=True)
    
    # Flatten multi-level columns for summary stats
    summary_flat = summary_stats.copy()
    summary_flat.columns = ['_'.join(col).strip() for col in summary_flat.columns.values]
    
    # Round for better readability
    summary_flat = summary_flat.round(4)
    
    # Save detailed results
    df.to_csv(f'{save_dir}/detailed_results.csv', index=False)
    summary_flat.to_csv(f'{save_dir}/summary_statistics.csv')
    
    # Create ranking table
    ranking_metrics = {
        'Best Final RMSD': 'final_best_fitness_mean',
        'Best Improvement': 'improvement_mean', 
        'Fastest Execution': 'execution_time_mean',
        'Most Consistent': 'final_best_fitness_std'
    }
    
    rankings = {}
    for metric_name, column in ranking_metrics.items():
        if column in summary_flat.columns:
            if 'time' in column or 'rmsd' in column.lower() or 'std' in column:
                # Lower is better
                rankings[metric_name] = summary_flat[column].rank().astype(int)
            else:
                # Higher is better
                rankings[metric_name] = summary_flat[column].rank(ascending=False).astype(int)
    
    ranking_df = pd.DataFrame(rankings, index=summary_flat.index)
    ranking_df['Average Rank'] = ranking_df.mean(axis=1).round(2)
    ranking_df = ranking_df.sort_values('Average Rank')
    
    ranking_df.to_csv(f'{save_dir}/configuration_rankings.csv')
    
    print("\n" + "="*60)
    print("CONFIGURATION RANKINGS")
    print("="*60)
    print(ranking_df.to_string())
    
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(summary_flat.to_string())
    
    return ranking_df, summary_flat

def run_comprehensive_test():
    """Run comprehensive EA testing"""
    print("Loading proteins...")
    
    # Create dummy protein data if files don't exist
    dummy_proteins = []
    for i in range(3):
        n_residues = np.random.randint(50, 150)
        coords = np.random.randn(n_residues, 3) * 10
        sequence = ''.join(np.random.choice(list('ACDEFGHIKLMNPQRSTVWY'), size=n_residues))
        
        dummy_proteins.append({
            'coords': coords,
            'sequence': sequence,
            'ca_atoms': None,
            'structure': None
        })
    
    proteins = dummy_proteins
    print(f"Loaded {len(proteins)} proteins")
    
    # Get test configurations
    configs = create_test_configs()
    print(f"Testing {len(configs)} configurations")
    
    # Run tests
    all_results = []
    all_histories = []
    
    for config in configs:
        results, histories = run_single_test(config, proteins, generations=30, trials=3)
        all_results.extend(results)
        all_histories.extend(histories)
    
    # Analyze results
    print("\nAnalyzing results...")
    df, summary_stats = analyze_results(all_results)
    
    # Create plots
    print("Creating comparison plots...")
    plot_comparison_results(all_results, all_histories)
    
    # Create tables
    print("Creating summary tables...")
    ranking_df, summary_flat = create_summary_tables(df, summary_stats)
    
    # Best configuration analysis
    best_config = ranking_df.index[0]
    print(f"\n{'='*60}")
    print(f"BEST OVERALL CONFIGURATION: {best_config}")
    print(f"{'='*60}")
    
    best_config_data = df[df['config'] == best_config]
    print(f"Mean Final RMSD: {best_config_data['final_best_fitness'].mean():.4f} ± {best_config_data['final_best_fitness'].std():.4f}")
    print(f"Mean Improvement: {best_config_data['improvement'].mean():.4f} ± {best_config_data['improvement'].std():.4f}")
    print(f"Mean Execution Time: {best_config_data['execution_time'].mean():.2f} ± {best_config_data['execution_time'].std():.2f} seconds")
    
    # Selection method analysis
    print(f"\n{'='*40}")
    print("SELECTION METHOD ANALYSIS")
    print(f"{'='*40}")
    
    # Extract selection method from config name more robustly
    def extract_selection_method(config_name):
        if 'Binary_Tournament' in config_name:
            return 'Binary_Tournament'
        elif 'Tournament' in config_name:
            return 'Tournament'
        elif 'Roulette' in config_name:
            return 'Roulette'
        elif 'Rank' in config_name:
            return 'Rank'
        else:
            return 'Unknown'
    
    df['selection_method'] = df['config'].apply(extract_selection_method)
    selection_analysis = df.groupby('selection_method')['final_best_fitness'].agg(['mean', 'std', 'count']).round(4)
    print(selection_analysis.to_string())
    
    print("\nTesting completed successfully!")
    print("Results saved in 'results/' directory")

if __name__ == "__main__":
    run_comprehensive_test()