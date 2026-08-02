# Conformation Dynamics: A Physics-Based Model of Protein State Propagation

A statistical physics and lattice model exploring non-genetic conformational state propagation in proteins, simulating induced misfolding mechanisms analogous to prion-like behavior.

![Propagation Dynamics](propagation_dynamics.png)

---

## 🔬 Core Scientific Question

How do physical contact interactions between protein conformations lower activation energy barriers to drive self-propagating misfolding across a spatial ensemble without genetic alterations?

---

## 📐 Mathematical & Physical Framework

### 1. Spontaneous Conformational Transition Rate
The probability of an isolated protein undergoing a spontaneous transition from Native state ($\text{State } C$) to Misfolded state ($\text{State } S$) per Monte Carlo step follows standard Boltzmann thermal activation:

$$P_{\text{spontaneous}} = \exp\left(-\frac{\Delta E}{k_B T}\right)$$

Where:
- $\Delta E$: Activation energy barrier ($\text{Joules}, \text{J}$)
- $k_B$: Boltzmann constant ($1.380649 \times 10^{-23} \text{ J/K}$)
- $T$: Absolute temperature ($\text{Kelvin}, \text{K}$)

---

### 2. Catalytic / Induced Transition Rate
Steric contact with neighboring misfolded proteins ($\text{State } S$) provides catalytic surface interaction, linearly reducing the effective activation energy barrier:

$$\Delta E_{\text{eff}} = \max\left(0, \Delta E - (\Delta E_{\text{lowering}} \cdot N_S)\right)$$

The catalytic transition rate constant is calculated as:

$$\text{Rate}_{\text{induced}} = k_{\text{cat}} \cdot N_S \cdot \exp\left(-\frac{\Delta E_{\text{eff}}}{k_B T}\right)$$

The resulting transition probability within a synchronous discrete-time step is:

$$P_{\text{induced}} = 1 - \exp\left(-\text{Rate}_{\text{induced}}\right)$$

Where:
- $N_S$: Number of nearest misfolded neighbors ($N_S \in \{0, 1, 2, 3, 4\}$ on a 2D square lattice).
- $\Delta E_{\text{lowering}}$: Barrier reduction constant per contact ($\text{J}$).
- $k_{\text{cat}}$: Catalytic rate scaling factor ($\text{step}^{-1}$).

---

### 3. Total Transition Probability
$$\text{Total Transition Probability } P_{\text{total}} = \min\left(1.0, P_{\text{spontaneous}} + P_{\text{induced}}\right)$$

---

## 🕹️ Simulation Topology & Update Rules
- **Lattice Dimensions:** $50 \times 50$ discrete 2D grid ($2,500$ total protein nodes).
- **Boundary Conditions:** 2D Periodic Boundary Conditions (Toroidal topology) to prevent edge reflection artifacts.
- **Update Scheme:** Synchronous state updating across all lattice sites per Monte Carlo Step.
- **Statistical Ensemble:** Multi-seed Monte Carlo simulations ($N=5$) to extract standard deviation bands ($\pm 1 \sigma$).

---

## 🚀 Quick Start & Installation

### Requirements
Ensure you have Python 3.8+ installed:
```bash
pip install -r requirements.txt
python main.py
.
├── core/
│   └── physics.py         # Thermodynamic transition probability equations
├── simulation/
│   └── grid.py            # 2D Lattice model & periodic boundary conditions
├── main.py                # Ensemble runner & visualization pipeline
├── propagation_dynamics.png # High-resolution output graph
├── requirements.txt       # Project dependencies
└── README.md              # Documentation
