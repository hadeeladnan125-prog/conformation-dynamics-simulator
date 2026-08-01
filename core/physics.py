import numpy as np

# Physical Constants
K_BOLTZMANN = 1.380649e-23  # J/K


def calculate_spontaneous_prob(delta_e, temperature):
    """Calculates spontaneous transition probability using Boltzmann distribution.

    Args:
        delta_e (float): Activation energy barrier (Joules).
        temperature (float): System temperature (Kelvin).

    Returns:
        float: Probability of spontaneous conversion (0 to 1).
    """
    if temperature <= 0:
        return 0.0
    kb_t = K_BOLTZMANN * temperature
    return np.exp(-delta_e / kb_t)


def calculate_induced_prob(
    delta_e, delta_e_lowering, k_cat, num_s_neighbors, temperature
):
    """Calculates induced (catalytic) transition probability.

    Args:
        delta_e (float): Original activation energy barrier (Joules).
        delta_e_lowering (float): Barrier reduction per misfolded neighbor
          (Joules).
        k_cat (float): Catalytic rate scaling factor.
        num_s_neighbors (int): Number of adjacent State S neighbors.
        temperature (float): System temperature (Kelvin).

    Returns:
        float: Probability of induced conversion (0 to 1).
    """
    if num_s_neighbors == 0 or temperature <= 0:
        return 0.0

    kb_t = K_BOLTZMANN * temperature
    effective_barrier = max(0.0, delta_e - (delta_e_lowering * num_s_neighbors))
    rate = k_cat * num_s_neighbors * np.exp(-effective_barrier / kb_t)

    # Convert rate to probability (bounded between 0 and 1)
    prob = 1.0 - np.exp(-rate)
    return float(np.clip(prob, 0.0, 1.0))
  
