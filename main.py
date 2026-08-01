import matplotlib.pyplot as plt
from simulation.grid import ProteinLattice

# Simulation Setup Parameters
GRID_SIZE = 50
STEPS = 100
TEMPERATURE = 310.15  # Physiological body temperature (37°C in Kelvin)
DELTA_E = 1.2e-19  # Base activation energy barrier (Joules)
DELTA_E_LOWERING = 0.3e-19  # Barrier lowering per misfolded neighbor
K_CAT = 1.5e18  # Catalytic scale factor


def run_simulation():
  print("=" * 60)
  print(
      "Starting Physics-Informed Conformation Dynamics Simulation..."
  )
  print(f"System Temp: {TEMPERATURE} K | Grid: {GRID_SIZE}x{GRID_SIZE}")
  print("=" * 60)

  lattice = ProteinLattice(size=GRID_SIZE, initial_prions=1)

  c_history = []
  s_history = []

  for step in range(STEPS):
    c_count, s_count = lattice.step(
        delta_e=DELTA_E,
        delta_e_lowering=DELTA_E_LOWERING,
        k_cat=K_CAT,
        temperature=TEMPERATURE,
    )
    c_history.append(c_count)
    s_history.append(s_count)

    if step % 10 == 0 or step == STEPS - 1:
      print(f"Step {step:03d} | State C (Normal): {c_count:4d} | State S"
            f" (Prion): {s_count:4d}")

  # Plotting results
  plt.figure(figsize=(8, 5))
  plt.plot(c_history, label="State C (Native)", color="blue")
  plt.plot(s_history, label="State S (Misfolded/Prion)", color="red")
  plt.title("Conformational State Propagation Dynamics")
  plt.xlabel("Monte Carlo Time Steps")
  plt.ylabel("Protein Population")
  plt.legend()
  plt.grid(True)
  plt.savefig("propagation_dynamics.png")
  print("\nSimulation complete! Graph saved as 'propagation_dynamics.png'.")


if __name__ == "__main__":
  run_simulation()
  
