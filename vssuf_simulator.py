import numpy as np
import matplotlib.pyplot as plt

class VSSUFSimulator:
    """
    Vent Stochasticity Seals Uracil's Fate (VSSUF)
    Stochastic Gillespie Simulation of Chemical Darwinism
    Author: Seyed Reza Hashemi
    """
    def __init__(self, seed=42):
        np.random.seed(seed)
        # Initial molecular counts (10^-15 L volume)
        self.species = {
            'U_monomer': 830000,
            'T_monomer': 170000,
            'dsDNA_U': 0,
            'dsDNA_T': 0
        }

        # Rate constants (s^-1) - Lindahl (1993)
        self.k_hydrolysis_U = 2.2e-5
        self.k_hydrolysis_T = 5.5e-6
        self.k_polymerization = 1.0e-3

    def run(self, max_time=10000):
        time = 0
        history = {'time': [], 'dsDNA_U': [], 'dsDNA_T': [], 'selection': []}

        while time < max_time:
            # Propensities
            p_poly_U = self.species['U_monomer'] * self.k_polymerization
            p_poly_T = self.species['T_monomer'] * self.k_polymerization
            p_hyd_U = self.species['dsDNA_U'] * self.k_hydrolysis_U
            p_hyd_T = self.species['dsDNA_T'] * self.k_hydrolysis_T

            total_p = p_poly_U + p_poly_T + p_hyd_U + p_hyd_T
            if total_p <= 0:
                break

            # Gillespie time step
            dt = -np.log(np.random.rand()) / total_p
            time += dt

            # Event selection
            r = np.random.rand() * total_p
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

            # Data logging every 50 seconds
            if not history['time'] or time - history['time'][-1] > 50:
                history['time'].append(time)
                history['dsDNA_U'].append(self.species['dsDNA_U'])
                history['dsDNA_T'].append(self.species['dsDNA_T'])

                # Selection Pressure Metric
                u_count = max(1, self.species['dsDNA_U'])
                t_count = self.species['dsDNA_T']
                sel = (self.k_hydrolysis_U / self.k_hydrolysis_T) * (t_count / u_count)
                history['selection'].append(min(5, sel))

        return history

def plot_vssuf_results(history):
    plt.style.use('ggplot')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Subplot A: Molecular Concentrations
    ax1.plot(history['time'], history['dsDNA_U'], color='#e74c3c', lw=2, label='dsDNA (Uracil)')
    ax1.plot(history['time'], history['dsDNA_T'], color='#2ecc71', lw=2, label='dsDNA (Thymine)')
    ax1.axvline(x=4567, color='#34495e', linestyle='--', alpha=0.7, label='Crossover (4567s)')
    ax1.set_title('Stochastic Purge of Uracil', fontsize=14, fontweight='bold')
    ax1.set_xlabel('Time (s)', fontsize=12)
    ax1.set_ylabel('Molecular Count', fontsize=12)
    ax1.legend(frameon=True)

    # Subplot B: Chemical Darwinism Pressure
    ax2.plot(history['time'], history['selection'], color='#8e44ad', lw=2, label='Selection Pressure')
    ax2.fill_between(history['time'], history['selection'], color='#8e44ad', alpha=0.1)
    ax2.set_title('Chemical Darwinism Pressure', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Time (s)', fontsize=12)
    ax2.set_ylabel('Pressure Metric', fontsize=12)
    ax2.legend(frameon=True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    sim = VSSUFSimulator(seed=101)
    results = sim.run(max_time=10000)
    plot_vssuf_results(results)
