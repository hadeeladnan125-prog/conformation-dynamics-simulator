import numpy as np
import matplotlib.pyplot as plt
from core.physics import set_random_seed, K_BOLTZMANN
from simulation.grid import ProteinLattice

def analytical_solution(steps, size, delta_e, temperature):
    """Calculates the theoretical baseline conversion assuming purely spontaneous transitions."""
    total_sites = size * size
    kb_t = K_BOLTZMANN * temperature
    p_spont = np.exp(-delta_e / kb_t)
    
    time_points = np.arange(steps + 1)
    analytical_s = total_sites * (1.0 - np.exp(-p_spont * time_points))
    return time_points, analytical_s

def run_multi_seed_simulation(size=50, steps=100, num_runs=5, delta_e=5e-21, 
                               delta_e_lowering=1e-21, k_cat=0.5, temp=310.15):
    """Runs multiple simulation passes to aggregate statistics (Mean and Std Dev)."""
    all_s_counts = []
    
    for seed in range(42, 42 + num_runs):
        set_random_seed(seed)
        lattice = ProteinLattice(size=size, initial_prions=1, seed=seed)
        
        s_history = [1]
        for _ in range(steps):
            _, s_cnt = lattice.step(delta_e, delta_e_lowering, k_cat, temp)
            s_history.append(s_cnt)
            
        all_s_counts.append(s_history)
        
    return np.array(all_s_counts)

def main():
    steps = 100
    grid_size = 50
    num_runs = 5
    temp = 310.15  # 37°C in Kelvin
    
    # Energy parameters in Joules
    delta_e = 5e-21
    delta_e_lowering = 1.2e-21
    k_cat = 0.6
    
    print("Running protein propagation ensemble simulations...")
    s_counts_matrix = run_multi_seed_simulation(
        size=grid_size, steps=steps, num_runs=num_runs, 
        delta_e=delta_e, delta_e_lowering=delta_e_lowering, 
        k_cat=k_cat, temp=temp
    )
    
    mean_s = np.mean(s_counts_matrix, axis=0)
    std_s = np.std(s_counts_matrix, axis=0)
    time_steps, analytical_s = analytical_solution(steps, grid_size, delta_e, temp)
    
    # Plotting Propagation Dynamics with Error Bands and Analytical Model
    fig, ax = plt.subplots(figsize=(9, 5), dpi=300)
    
    ax.plot(time_steps, mean_s, color='#d62728', lw=2.5, label='Lattice Model (Mean)')
    ax.fill_between(time_steps, mean_s - std_s, mean_s + std_s, color='#d62728', alpha=0.2, label='±1 Std Dev (Ensemble)')
    ax.plot(time_steps, analytical_s, color='#1f77b4', linestyle='--', lw=2.0, label='Purely Spontaneous Baseline (Analytical)')
    
    ax.set_title("Non-Genetic Conformational State Propagation in Proteins", fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel("Monte Carlo Steps (Time)", fontsize=11)
    ax.set_ylabel("Misfolded Population (State S Count)", fontsize=11)
    ax.set_xlim(0, steps)
    ax.set_ylim(0, grid_size * grid_size)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.legend(frameon=True, loc='upper left')
    
    plt.tight_layout()
    plt.savefig("propagation_dynamics.png", dpi=300)
    print("Simulation complete. Saved high-res figure 'propagation_dynamics.png'.")

if __name__ == "__main__":
    main()
  
