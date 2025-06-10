import numpy as np
from features import FEATURE_FUNCTIONS
from utils import rmsd
from tqdm import tqdm

class ParticleSwarmOptimizer:
    def _init_(self, proteins, population_size=30, generations=50, inertia=0.7, cognitive=1.5, social=1.5):
        self.proteins = proteins
        self.population_size = population_size
        self.generations = generations
        self.inertia = inertia
        self.cognitive = cognitive
        self.social = social
        self.vector_size = len(FEATURE_FUNCTIONS)

        # Particle positions and velocities
        self.positions = [self._init_particle() for _ in range(self.population_size)]
        self.velocities = [np.random.randn(self.vector_size) * 0.1 for _ in range(self.population_size)]
        self.personal_best = self.positions[:]
        self.personal_best_scores = [self.evaluate(p) for p in self.positions]
        self.global_best = self.personal_best[np.argmin(self.personal_best_scores)]

    def _init_particle(self):
        return np.random.randn(self.vector_size)

    def vector_to_coords(self, vector, target_length):
        padded = np.tile(vector, (target_length, 1))[:, :3]
        return padded

    def evaluate(self, weights):
        total_rmsd = 0
        for protein in self.proteins:
            features = [f(protein) for f in FEATURE_FUNCTIONS]
            weighted = sum(w * f for w, f in zip(weights, features))
            pred_coords = self.vector_to_coords(weighted, len(protein['coords']))
            total_rmsd += rmsd(pred_coords, protein['coords'])
        return total_rmsd / len(self.proteins)

    def run(self):
        for gen in range(self.generations):
            for i in range(self.population_size):
                r1, r2 = np.random.rand(), np.random.rand()
                cognitive_term = self.cognitive * r1 * (self.personal_best[i] - self.positions[i])
                social_term = self.social * r2 * (self.global_best - self.positions[i])
                self.velocities[i] = self.inertia * self.velocities[i] + cognitive_term + social_term
                self.positions[i] += self.velocities[i]

                score = self.evaluate(self.positions[i])
                if score < self.personal_best_scores[i]:
                    self.personal_best_scores[i] = score
                    self.personal_best[i] = self.positions[i].copy()

            self.global_best = self.personal_best[np.argmin(self.personal_best_scores)]
            best_score = min(self.personal_best_scores)
            print(f"Generation {gen}, Best RMSD: {best_score:.3f}")
            