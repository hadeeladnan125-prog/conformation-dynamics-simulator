import numpy as np

# Physical Constants
K_BOLTZMANN = 1.380649e-23  # Joules per Kelvin (J/K)


def set_random_seed(seed=42):
    """Fixes the random seed for reproducibility across simulation runs."""
    np.random.seed(seed)


def calculate_spontaneous_prob(delta_e, temperature):
    """Calculates the spontaneous conformational transition probability.

    Uses the standard Boltzmann factor based on thermal activation energy.

    Equation:
        P_spontaneous = exp(-Delta_E / (k_B * T))

    Args:
        delta_e (float): Activation energy barrier in Joules (J).
        temperature (float): Absolute temperature in Kelvin (K).

    Returns:
        float: Spontaneous transition probability within [0, 1].
    """
    if temperature <= 0:
        return 0.0
    kb_t = K_BOLTZMANN * temperature
    return float(np.exp(-delta_e / kb_t))


def calculate_induced_prob(
    delta_e, delta_e_lowering, k_cat, num_s_neighbors, temperature
):
    """Calculates the catalytic (induced) transition probability.

    Steric contact with misfolded neighbors reduces the effective activation
    energy barrier linearly.

    Equation:
        Effective_Barrier = max(0, Delta_E - (delta_e_lowering * N_S))
        Rate = k_cat * N_S * exp(-Effective_Barrier / (k_B * T))
        P_induced = 1 - exp(-Rate)

    Args:
        delta_e (float): Uncatalyzed activation energy barrier in Joules (J).
        delta_e_lowering (float): Barrier reduction per misfolded neighbor in
          Joules (J).
        k_cat (float): Catalytic scaling rate constant (1/step).
        num_s_neighbors (int): Number of adjacent State S (misfolded) neighbors
          (0 to 4).
        temperature (float): Absolute temperature in Kelvin (K).

    Returns:
        float: Induced conversion probability within [0, 1].
    """
    if num_s_neighbors == 0 or temperature <= 0:
        return 0.0

    kb_t = K_BOLTZMANN * temperature
    effective_barrier = max(0.0, delta_e - (delta_e_lowering * num_s_neighbors))
    rate = k_cat * num_s_neighbors * np.exp(-effective_barrier / kb_t)

    prob = 1.0 - np.exp(-rate)
    return float(np.clip(prob, 0.0, 1.0))
    
