import numpy as np
from core.physics import calculate_induced_prob, calculate_spontaneous_prob

# State Definitions
STATE_C = 1  # Normal Protein
STATE_S = 2  # Misfolded / Prion Protein


class ProteinLattice:

    def __init__(self, size=50, initial_prions=1):
        """Initializes the 2D Protein Lattice system."""
        self.size = size
        self.grid = np.ones((size, size), dtype=int) * STATE_C

        # Plant initial prion seeds at the center
        center = size // 2
        for i in range(initial_prions):
            self.grid[center + i, center] = STATE_S

    def count_s_neighbors(self, r, c):
        """Counts State S neighbors with periodic boundary conditions."""
        neighbors = [
            self.grid[(r - 1) % self.size, c],
            self.grid[(r + 1) % self.size, c],
            self.grid[r, (c - 1) % self.size],
            self.grid[r, (c + 1) % self.size],
        ]
        return neighbors.count(STATE_S)

    def step(self, delta_e, delta_e_lowering, k_cat, temperature):
        """Advances the simulation by one Monte Carlo time step."""
        new_grid = self.grid.copy()

        for r in range(self.size):
            for c in range(self.size):
                if self.grid[r, c] == STATE_C:
                    s_neighbors = self.count_s_neighbors(r, c)

                    p_spont = calculate_spontaneous_prob(delta_e, temperature)
                    p_ind = calculate_induced_prob(
                        delta_e, delta_e_lowering, k_cat, s_neighbors, temperature
                    )

                    # Total probability of conversion
                    p_total = np.clip(p_spont + p_ind, 0.0, 1.0)

                    if np.random.rand() < p_total:
                        new_grid[r, c] = STATE_S

        self.grid = new_grid
        c_count = np.sum(self.grid == STATE_C)
        s_count = np.sum(self.grid == STATE_S)
        return c_count, s_count
      
