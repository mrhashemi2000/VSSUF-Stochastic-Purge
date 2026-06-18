
# Vent Stochasticity Seals Uracil's Fate: A Physical Model of Chemical Darwinism in the Origin of DNA

Author: Seyed Mohammad Reza Hashemi (Reza Hashemi)
ORCID: 0009-0002-0645-5180  
Affiliation: Former Member of the Pasteur Institute of Iran  
Repository: VSSUF-Stochastic-Purge
Abstract
This project demonstrates that the exclusion of uracil from DNA in favor of thymine is an inevitable outcome of prebiotic physical chemistry. Using stochastic Gillespie simulations, we model the fluctuating conditions of hydrothermal vents—including diffusion-limited encounters, Poisson-distributed ROS bursts, and discrete molecular interactions. 

Key Result: Over a simulated vent lifetime, the fraction of uracil in double-stranded DNA collapses from 83% to 7.2% (n=50 runs; p=0.002). This chemical selection is driven by thymine's superior hydrolytic stability and replication fidelity.Scientific Core: Chemical Darwinism
This model provides physical validation for "Chemical Darwinism": the concept that Darwinian selection operates at the molecular level as a direct consequence of physicochemical laws. The system evolves through four distinct phases:
1. Promiscuous Polymerization: Dominance of uracil.
2. Duplex Emergence: Thymine's structural compatibility leads to faster annealing (2.4x).
3. The Purge: Rapid elimination of uracil once selection pressure exceeds a critical threshold (>2.0).
4. Canonical Lock-in: Stabilization of the thymine-rich genome.Methodology & Implementation
The simulation is built on a rigorous physical foundation:
- Algorithm: Custom implementation of the Gillespie Stochastic Simulation Algorithm (SSA).
- Empirical Data: Hydrolysis rates are derived from Lindahl (1993).
- Environmental Parameters: Base temperature 80^C at pH 6.2.
- Robustness: The "Uracil Purge" is robust to pm 20 changes in ROS flux and 50% reductions in diffusion rates.
  How to Run the Simulation
1. Clone the repository:
   git clone https://github.com/mrhashemi2000/VSSUF-Stochastic-Purge.git
   2. Install dependencies:
   pip install -r requirements.txt
Run the simulator:
  bash
   python vssuf_simulator.py
References
- Lindahl, T. (1993). *Nature* 362, 709–715.
- Ranjan, S. et al. (2019). *Nat. Astron.* 3, 140–142.
- Becker, S. et al. (2019). *Science* 366, 76–82.
- Hashemi, R. (2025). The Matter World Hypothesis: Chemical Darwinism Computational Simulation Origin of Life I-VII. Zenodo. https://doi.org/10.5281/zenodo.17650234
DOI: https://doi.org/10.5281/zenodo.17273763
