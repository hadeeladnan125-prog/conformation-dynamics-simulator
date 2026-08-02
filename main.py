import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

from core.physics import K_BOLTZMANN, set_random_seed
from simulation.grid import ProteinLattice

GRID_SIZE = 50
TOTAL_NODES = GRID_SIZE * GRID_SIZE
STEPS = 100
RUNS = 20
TEMPERATURE = 310.15  # 37°C in Kelvin
DELTA_E = 1.2e-19  # Joules
DELTA_E_LOWERING = 0.3e-19  # Joules
K_CAT = 1.5e18


def logistic_model(t, k, t0):
  """Analytical Mean-Field Logistic Differential Equation Solution:

  S(t) = N_total / (1 + exp(-k * (t - t0)))
  """
  return TOTAL_NODES / (1.0 + np.exp(-k * (t - t0)))


def run_statistical_simulation():
  print("=" * 60)
  print(f"Running Statistical Ensemble ({RUNS} runs) & Validation...")
  print("=" * 60)

  all_s_curves = np.zeros((RUNS, STEPS + 1))

  for run in range(RUNS):
    seed = 42 + run
    set_random_seed(seed)
    lattice = ProteinLattice(size=GRID_SIZE, initial_prions=1, seed=seed)

    s_history = [1]
    for step in range(STEPS):
      _, s_count = lattice.step(
          DELTA_E, DELTA_E_LOWERING, K_CAT, TEMPERATURE
      )
      s_history.append(s_count)

    all_s_curves[run, :] = s_history

  # Statistical Aggregation
  s_mean = np.mean(all_s_curves, axis=0)
  s_std = np.std(all_s_curves, axis=0)
  time_steps = np.arange(STEPS + 1)

  # Analytical Mean-Field Curve Fit
  popt, _ = curve_fit(
      logistic_model, time_steps, s_mean, p0=[0.1, STEPS / 2]
  )
  s_analytical = logistic_model(time_steps, *popt)

  # Quantitative Validation Metric (R-squared)
  ss_res = np.sum((s_mean - s_analytical) ** 2)
  ss_tot = np.sum((s_mean - np.mean(s_mean)) ** 2)
  r_squared = 1 - (ss_res / ss_tot)

  # Plotting Simulation vs Analytical Fit
  fig, ax = plt.subplots(figsize=(9, 5.5), dpi=300)

  # Ensemble Simulation Curve
  ax.plot(
      time_steps,
      s_mean,
      color="#d62728",
      linewidth=2.5,
      label=f"2D Lattice Simulation (Mean of {RUNS} runs)",
  )
  ax.fill_between(
      time_steps,
      s_mean - s_std,
      s_mean + s_std,
      color="#d62728",
      alpha=0.25,
      label="±1 Std Deviation (Ensemble)",
  )

  # Analytical Fit
  ax.plot(
      time_steps,
      s_analytical,
      color="#002b36",
      linestyle="--",
      linewidth=2.0,
      label=f"Analytical Mean-Field Fit ($R^2 = {r_squared:.4f}$)",
  )

  ax.set_title(
      "Spatial Lattice Simulation vs. Analytical Mean-Field Model",
      fontsize=12,
      fontweight="bold",
      pad=12,
  )
  ax.set_xlabel("Monte Carlo Time Steps", fontsize=11)
  ax.set_ylabel("Misfolded Population (State S Count)", fontsize=11)
  ax.set_xlim(0, STEPS)
  ax.set_ylim(0, TOTAL_NODES)
  ax.grid(True, linestyle=":", alpha=0.6)
  ax.legend(loc="upper left", frameon=True)

  plt.tight_layout()
  plt.savefig("propagation_dynamics.png", dpi=300)
  print(
      f"Simulation complete! Graph saved as 'propagation_dynamics.png' (R^2 ="
      f" {r_squared:.4f})."
  )


if __name__ == "__main__":
  run_statistical_simulation()
  
