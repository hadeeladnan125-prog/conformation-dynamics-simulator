# Conformation Dynamics: A Physics-Based Model of Protein State Propagation

## Core Scientific Question
> **Can structural conformation transition alone, governed by fundamental physical laws, generate self-propagating infection-like dynamics without the presence of DNA or RNA?**

This project presents a physics-informed lattice simulation model to study non-genetic conformational state propagation in proteins (prion dynamics).

---

## Model Architecture

The simulation treats protein units as lattice nodes operating under a two-state conformational paradigm:

- **State C (Normal Conformation):** Thermodynamically metastable state representing native protein structure.
- **State S (Mis-folded / Prion Conformation):** Globally stable state with lower free energy ($E_S < E_C$), but separated by an activation energy barrier $\Delta E$.
### Physical Parameters & Governing Laws
1. **Boltzmann Transition Probability:**
   $$P_{\text{spontaneous}} = \exp\left(-\frac{\Delta E}{k_B T}\right)$$
2. **Catalytic Transition (Induced Conversion):**
   $$P_{\text{induced}} = 1 - \exp\left(-k_{\text{cat}} \cdot N_S \cdot \exp\left(-\frac{\Delta E - \delta E}{k_B T}\right)\right)$$
   where $N_S$ is the number of adjacent misfolded neighbors, and $\delta E$ represents the activation barrier lowering due to steric contact.
### Simulation Results
Below is the numerical result showing the sigmoidal propagation of the misfolded state $S$ over time:

![Conformation Propagation Dynamics](propagation_dynamics.png)
---

## Project Structure

```text
conformation-dynamics-simulator/
├── core/
│   ├── __init__.py
│   └── physics.py         # Thermodynamic transition rate functions
├── simulation/
│   ├── __init__.py
│   └── grid.py            # 2D Lattice simulation engine & state tracking
├── main.py                # Main simulation entry point & execution script
├── requirements.txt       # Python dependencies
└── README.md              # Project documentation
git clone [https://github.com/hadeeladnan125-prog/conformation-dynamics-simulator.git](https://github.com/hadeeladnan125-prog/conformation-dynamics-simulator.git)
cd conformation-dynamics-simulator
pip install -r requirements.txt
python main.py

