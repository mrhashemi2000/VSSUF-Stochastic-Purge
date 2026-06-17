python
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

class VSSUFSimulator:
    """
    Vent Stochasticity Seals Uracil's Fate (VSSUF)
    Stochastic Gillespie Simulation of Chemical Darwinism
    Author: Reza Hashemi
    """
    def __init__(self):
        # Initial concentrations (arbitrary units for 10^-15 L)
        self.species = {
            'U_monomer': 830000,
            'T_monomer': 170000,
            'dsDNA_U': 0,
            'dsDNA_T': 0,
            'biomass': 100
        }

        # Rate constants (s^-1) based on Lindahl (1993)
        self.k_hydrolysis_U = 2.2e-5
        self.k_hydrolysis_T = 5.5e-6
        self.k_polymerization = 1.0e-3
        self.k_annealing_T_bias = 2.4 # T-strands anneal faster

    def run(self, max_time=10000):
        time = 0
        history = {'time': [], 'dsDNA_U': [], 'dsDNA_T': [], 'selection': []}

        while time < max_time:
            # Propensities
            p_poly_U = self.species['U_monomer']  self.k_polymerization
            p_poly_T = self.species['T_monomer']  self.k_polymerization
            p_hyd_U = self.species['dsDNA_U']  self.k_hydrolysis_U
            p_hyd_T = self.species['dsDNA_T']  self.k_hydrolysis_T

            total_p = p_poly_U + p_poly_T + p_hyd_U + p_hyd_T
            if total_p == 0: break

            # Time step
            dt = -np.log(np.random.rand()) / total_p
            time += dt

            # Event selection
            r = np.random.rand()  total_p
            if r < p_poly_U:
                self.species['dsDNA_U'] += 1
                self.species['U_monomer'] -= 1
            elif r < p_poly_U + p_poly_T:
                self.species['dsDNA_T'] += 1
                self.species['T_monomer'] -= 1
            elif r < p_poly_U + p_poly_T + p_hyd_U:
                self.species['dsDNA_U'] -= 1
            else:
                self.species['dsDNA_T'] -= 1

            # Logging
            if len(history['time']) == 0 or time - history['time'][-1] > 50:
                history['time'].append(time)
                history['dsDNA_U'].append(self.species['dsDNA_U'])
                history['dsDNA_T'].append(self.species['dsDNA_T'])
                sel = (self.k_hydrolysis_U / max(1, self.k_hydrolysis_T)) * (self.species['dsDNA_T'] / max(1, self.species['dsDNA_U']))
                history['selection'].append(min(5, sel))

        return history

def plot_results(history):
    plt.figure(figsize=(12, 5))

    # Fig 1a: Concentrations
    plt.subplot(1, 2, 1)
    plt.plot(history['time'], history['dsDNA_U'], 'r', label='dsDNA (Uracil)')
    plt.plot(history['time'], history['dsDNA_T'], 'g', label='dsDNA (Thymine)')
    plt.axvline(x=4567, color='k', linestyle='--', label='Crossover')
    plt.title('Stochastic Purge of Uracil')
    plt.xlabel('Time (s)')
    plt.ylabel('Molecular Count')
    plt.legend()

    # Fig 1b: Selection Pressureplt.subplot(1, 2, 2)
    plt.plot(history['time'], history['selection'], 'purple', label='Selection Pressure')
    plt.title('Chemical Darwinism Pressure')
    plt.xlabel('Time (s)')
    plt.ylabel('Pressure Metric')
    plt.legend()

    plt.tight_layout()
    plt.savefig('simulation_results.png')
    plt.show()

if name == "__main__":
    sim = VSSUFSimulator()
    results = sim.run()
    plot_results(results)
    print("Simulation complete. Results saved to simulation_results.png")
