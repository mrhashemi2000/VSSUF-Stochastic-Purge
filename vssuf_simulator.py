import numpy as np
import matplotlib.pyplot as plt

class VSSUFSimulator:
    """
    Vent Stochasticity Seals Uracil's Fate (VSSUF)
    Physical Model: Arrhenius-calibrated Gillespie Simulation
    Author: Seyed Mohammad Reza Hashemi (Reza Hashemi)
    """
    def __init__(self, seed=101, temperature_C=80):
        # Use a fixed seed for reproducibility
        np.random.seed(seed)

        # --- Arrhenius Calibration Section ---
        self.T_kelvin = temperature_C + 273.15
        R = 8.314            # Gas constant
        T_ref = 310.15       # Reference temperature (37C)
        Ea = 105000.0        # Activation energy for hydrolysis

        # Base rates at 37C
        k_U_ref = 3.0e-8 
        k_T_ref = 7.5e-9 

        # Calculate Arrhenius factor for transfer to vent temperature
        inv_T_diff = (1.0 / T_ref) - (1.0 / self.T_kelvin)
        arrhenius_factor = np.exp((Ea / R) * inv_T_diff)

        self.k_hydrolysis_U = k_U_ref * arrhenius_factor
        self.k_hydrolysis_T = k_T_ref * arrhenius_factor

        print(f"--- VSSUF Calibration at {temperature_C}C ---")
        print(f"k_U: {self.k_hydrolysis_U:.2e} s^-1")
        print(f"k_T: {self.k_hydrolysis_T:.2e} s^-1")

        # --- Initial Molecular Conditions ---
        self.species = {
            'U_monomer': 830000, 
            'T_monomer': 170000,
            'dsDNA_U': 0,
            'dsDNA_T': 0
        }
        self.k_polymerization = 0.001

    def run(self, max_time=10000):
        time = 0.0
        history = {'time': [], 'dsDNA_U': [], 'dsDNA_T': [], 'selection': []}

        while time < max_time:
            # Calculate Gillespie Propensities (Reaction Probabilities)
            p_poly_U = self.species['U_monomer'] * self.k_polymerization
            p_poly_T = self.species['T_monomer'] * self.k_polymerization
            p_hyd_U = self.species['dsDNA_U'] * self.k_hydrolysis_U
            p_hyd_T = self.species['dsDNA_T'] * self.k_hydrolysis_T

            total_p = p_poly_U + p_poly_T + p_hyd_U + p_hyd_T
            if total_p <= 0:
                break

            # Stochastic Time Step (Determining time to next event)
            dt = -np.log(np.random.rand()) / total_p
            time += dt

            # Reaction Selection (Choosing which event occurs)
            r = np.random.rand() * total_p
            if r < p_poly_U:
                self.species['dsDNA_U'] += 1
                self.species['U_monomer'] -= 1
            elif r < (p_poly_U + p_poly_T):
                self.species['dsDNA_T'] += 1
                self.species['T_monomer'] -= 1
            elif r < (p_poly_U + p_poly_T + p_hyd_U):
                self.species['dsDNA_U'] -= 1
            else:
                self.species['dsDNA_T'] -= 1

            # Data Logging (Sample every 40 seconds)
            if not history['time'] or (time - history['time'][-1] > 40):
                history['time'].append(time)
                history['dsDNA_U'].append(self.species['dsDNA_U'])
                history['dsDNA_T'].append(self.species['dsDNA_T'])

                # Chemical Darwinism: Selection Pressure Calculation
                u_eff = max(1, self.species['dsDNA_U'])
                t_eff = self.species['dsDNA_T']
                sel_val = (self.k_hydrolysis_U / self.k_hydrolysis_T) * (t_eff / (u_eff + t_eff + 1e-9)) * 5.0
                history['selection'].append(min(5.0, sel_val))

        return history

def plot_results(history):
    plt.style.use('ggplot')
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Figure 1a: Molecular Purge of Uracil
    ax1.plot(history['time'], history['dsDNA_U'], color='#e74c3c', lw=2, label='dsDNA (Uracil)')
    ax1.plot(history['time'], history['dsDNA_T'], color='#2ecc71', lw=2, label='dsDNA (Thymine)')
    ax1.set_title('Figure 1a: Molecular Purge of Uracil', fontweight='bold')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Molecular Count')
    ax1.legend()

    # Figure 1b: Chemical Darwinism Pressure
    ax2.plot(history['time'], history['selection'], color='#8e44ad', lw=2, label='Selection Pressure')
    ax2.fill_between(history['time'], history['selection'], color='#8e44ad', alpha=0.1)
    ax2.set_title('Figure 1b: Chemical Darwinism Pressure', fontweight='bold')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Pressure Metric (Sp)')
    ax2.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    print("Initializing VSSUF Simulation (Gillespie SSA)...")
    sim = VSSUFSimulator(seed=101, temperature_C=80)
    results = sim.run(max_time=10000)
    plot_results(results)
