"""
Integrated Formate Beads + Monod Kinetics Calculator

This module combines formate bead release kinetics with bacterial consumption
using Monod kinetics to predict substrate concentration over time.
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from formate_beads_core import ExperimentManager, M07_BEAD_RELEASE, M03_BEAD_RELEASE, Bead
from monod_kinetics import MonodKinetics


class ConstantSubstrateCalculator:
    """
    Calculates bead schedule to maintain constant substrate concentration
    despite bacterial consumption.
    """
    
    def __init__(self, volume_ml, mu_max, Ks, Y_xs, target_substrate_mM, od_conversion=0.4):
        """
        Initialize calculator for constant substrate maintenance.
        
        Parameters:
        -----------
        volume_ml : float
            Culture volume (mL)
        mu_max : float
            Maximum specific growth rate (1/hr)
        Ks : float
            Half-saturation constant (mM)
        Y_xs : float
            Yield coefficient (g biomass / mmol substrate)
        target_substrate_mM : float
            Target substrate concentration to maintain (mM)
        od_conversion : float
            OD to biomass conversion (g/L per OD unit)
        """
        self.volume_ml = volume_ml
        self.volume_L = volume_ml / 1000.0
        self.monod = MonodKinetics(mu_max, Ks, Y_xs, volume_ml)
        self.target_substrate_mM = target_substrate_mM
        self.od_conversion = od_conversion
        self.beads = []
        self.bead_schedule = {}
    
    def calculate_consumption_rate(self, od):
        """
        Calculate bacterial substrate consumption rate at target concentration.
        
        Parameters:
        -----------
        od : float
            Current optical density
        
        Returns:
        --------
        float : Consumption rate (mmol/day)
        """
        # Convert OD to biomass
        biomass_g_L = self.monod.od_to_biomass(od, self.od_conversion)
        
        # Get growth rate at target substrate concentration
        mu = self.monod.specific_growth_rate(self.target_substrate_mM)
        
        # Calculate substrate uptake rate (mmol/(L·hr))
        uptake_rate_per_hr = (mu * biomass_g_L) / self.monod.Y_xs
        
        # Convert to mmol/day for the entire volume
        uptake_mmol_per_day = uptake_rate_per_hr * 24.0 * self.volume_L
        
        return uptake_mmol_per_day
    
    def calculate_bead_schedule(self, initial_od, od_trajectory, days=7):
        """
        Calculate bead schedule to maintain constant substrate concentration.
        
        Parameters:
        -----------
        initial_od : float
            Initial optical density
        od_trajectory : dict or callable
            Either dict {day: expected_OD} or function od_trajectory(day)
            Predicts bacterial density over time
        days : int
            Number of days to calculate
        
        Returns:
        --------
        dict : Bead schedule {day: {'M07': count, 'M03': count}}
        """
        print("\n" + "="*80)
        print("CONSTANT SUBSTRATE CONCENTRATION - BEAD CALCULATOR")
        print("="*80)
        print(f"\nTarget substrate concentration: {self.target_substrate_mM:.2f} mM")
        print(f"Culture volume: {self.volume_ml} mL ({self.volume_L:.3f} L)")
        print(f"Bacterial parameters:")
        print(f"  μmax = {self.monod.mu_max:.3f} hr⁻¹")
        print(f"  Ks = {self.monod.Ks:.2f} mM")
        print(f"  Y_x/s = {self.monod.Y_xs:.3f} g/mmol")
        print("\n" + "="*80 + "\n")
        
        self.bead_schedule = {}
        self.beads = []
        
        for day in range(1, days + 1):
            # Get expected OD for this day
            if callable(od_trajectory):
                expected_od = od_trajectory(day)
            else:
                expected_od = od_trajectory.get(day, initial_od)
            
            # Calculate bacterial consumption rate
            consumption_mmol_per_day = self.calculate_consumption_rate(expected_od)
            
            # Calculate release from existing beads
            existing_release_mmol_per_day = 0
            for bead in self.beads:
                release = bead.get_release_rate(day)
                existing_release_mmol_per_day += release
            
            # Additional release needed
            additional_needed_mmol_per_day = consumption_mmol_per_day - existing_release_mmol_per_day
            
            print(f"Day {day}:")
            print(f"  Expected OD:              {expected_od:.3f}")
            print(f"  Bacterial consumption:    {consumption_mmol_per_day:.3f} mmol/day")
            print(f"  Existing bead release:    {existing_release_mmol_per_day:.3f} mmol/day")
            print(f"  Additional needed:        {additional_needed_mmol_per_day:.3f} mmol/day")
            
            m07_count = 0
            m03_count = 0
            
            if additional_needed_mmol_per_day > 0:
                # Use M03 beads for sustained release (Day 1 release rate)
                m03_day1_release = M03_BEAD_RELEASE[1]
                m07_day1_release = M07_BEAD_RELEASE[1]
                
                # Strategy: use M03 for baseline, M07 for extra boost if needed
                m03_beads_float = additional_needed_mmol_per_day / m03_day1_release
                
                if m03_beads_float > 3:
                    # Use 2 M03 beads and fill gap with M07
                    m03_count = 2
                    remaining = additional_needed_mmol_per_day - (m03_count * m03_day1_release)
                    m07_count = max(0, int(round(remaining / m07_day1_release)))
                else:
                    # Just use M03 beads
                    m03_count = max(0, int(round(m03_beads_float)))
                    m07_count = 0
                
                # Calculate actual release
                actual_additional = (m07_count * m07_day1_release + 
                                   m03_count * m03_day1_release)
                actual_total_release = existing_release_mmol_per_day + actual_additional
                
                # Net balance
                net_balance = actual_total_release - consumption_mmol_per_day
                
                print(f"  → ADD {m07_count} M07 beads + {m03_count} M03 beads")
                print(f"  Total release:            {actual_total_release:.3f} mmol/day")
                print(f"  Net substrate balance:    {net_balance:+.3f} mmol/day")
                
                if abs(net_balance) < 0.1:
                    print(f"  ✓ Balanced (substrate stable)")
                elif net_balance > 0:
                    print(f"  ⚠ Slight accumulation ({net_balance:.3f} mmol/day)")
                else:
                    print(f"  ⚠ Slight depletion ({net_balance:.3f} mmol/day)")
                
                # Add new beads
                for _ in range(m07_count):
                    self.beads.append(Bead(day, 'M07', M07_BEAD_RELEASE))
                for _ in range(m03_count):
                    self.beads.append(Bead(day, 'M03', M03_BEAD_RELEASE))
            else:
                print(f"  → No beads needed (existing beads sufficient)")
                print(f"  Excess release:           {abs(additional_needed_mmol_per_day):.3f} mmol/day")
            
            self.bead_schedule[day] = {'M07': m07_count, 'M03': m03_count}
            print()
        
        # Summary
        total_m07 = sum(day['M07'] for day in self.bead_schedule.values())
        total_m03 = sum(day['M03'] for day in self.bead_schedule.values())
        
        print("="*80)
        print("SUMMARY")
        print("="*80)
        print(f"Total M07 beads: {total_m07}")
        print(f"Total M03 beads: {total_m03}")
        print(f"Total beads:     {total_m07 + total_m03}")
        print(f"\nThis schedule should maintain ~{self.target_substrate_mM:.2f} mM substrate")
        print("despite bacterial consumption over {days} days.")
        print("="*80 + "\n")
        
        return self.bead_schedule
    
    def estimate_od_trajectory(self, initial_od, time_days):
        """
        Estimate OD trajectory assuming constant substrate concentration.
        
        Parameters:
        -----------
        initial_od : float
            Initial optical density
        time_days : int
            Number of days
        
        Returns:
        --------
        dict : {day: estimated_OD}
        """
        # At constant substrate, growth rate is constant
        mu = self.monod.specific_growth_rate(self.target_substrate_mM)
        
        # Calculate doubling time
        if mu > 0:
            doubling_time_hr = np.log(2) / mu
        else:
            doubling_time_hr = np.inf
        
        print(f"\nEstimated growth at {self.target_substrate_mM:.2f} mM substrate:")
        print(f"  Growth rate: {mu:.3f} hr⁻¹")
        print(f"  Doubling time: {doubling_time_hr:.1f} hours")
        
        trajectory = {}
        for day in range(1, time_days + 1):
            time_hr = day * 24
            # Exponential growth: OD(t) = OD0 * e^(μ*t)
            od = initial_od * np.exp(mu * time_hr)
            trajectory[day] = od
        
        return trajectory


class IntegratedBeadBacteriaModel:
    """
    Models both formate release from beads AND bacterial consumption.
    """
    
    def __init__(self, volume_ml, mu_max, Ks, Y_xs, od_conversion=0.4):
        """
        Initialize integrated model.
        
        Parameters:
        -----------
        volume_ml : float
            Culture volume (mL)
        mu_max : float
            Maximum specific growth rate (1/hr)
        Ks : float
            Half-saturation constant (mM)
        Y_xs : float
            Yield coefficient (g biomass / mmol substrate)
        od_conversion : float
            OD to biomass conversion (g/L per OD unit)
        """
        self.volume_ml = volume_ml
        self.volume_L = volume_ml / 1000.0
        self.monod = MonodKinetics(mu_max, Ks, Y_xs, volume_ml)
        self.od_conversion = od_conversion
        self.bead_schedule = None
        self.beads = []
    
    def set_bead_schedule(self, schedule_dict):
        """
        Set when beads are added to the culture.
        
        Parameters:
        -----------
        schedule_dict : dict
            Dictionary mapping day to {'M07': count, 'M03': count}
            Example: {1: {'M07': 2, 'M03': 1}, 2: {'M07': 0, 'M03': 1}}
        """
        self.bead_schedule = schedule_dict
        
        # Create bead objects
        from formate_beads_core import Bead
        self.beads = []
        
        for day, counts in schedule_dict.items():
            for _ in range(counts['M07']):
                self.beads.append(Bead(day, 'M07', M07_BEAD_RELEASE))
            for _ in range(counts['M03']):
                self.beads.append(Bead(day, 'M03', M03_BEAD_RELEASE))
    
    def get_bead_release_rate(self, time_hours):
        """
        Calculate total formate release rate from all beads at given time.
        
        Parameters:
        -----------
        time_hours : float
            Time in hours
        
        Returns:
        --------
        float : Release rate (mmol/(L·hr))
        """
        if not self.beads:
            return 0.0
        
        # Convert hours to days
        current_day = int(time_hours / 24) + 1
        
        # Get total release for this day (mmol/day per bead)
        total_release_mmol_per_day = 0
        for bead in self.beads:
            release = bead.get_release_rate(current_day)
            total_release_mmol_per_day += release
        
        # Convert to mmol/(L·hr)
        # Release is mmol/day total, need to divide by volume and convert to per hour
        release_rate = total_release_mmol_per_day / (24.0 * self.volume_L)
        
        return release_rate
    
    def integrated_model(self, state, t):
        """
        Combined differential equations for bead release + bacterial consumption.
        
        Parameters:
        -----------
        state : tuple
            (biomass_g_L, substrate_mM)
        t : float
            Time (hours)
        
        Returns:
        --------
        tuple : (dX/dt, dS/dt)
        """
        biomass, substrate = state
        
        # Ensure non-negative
        biomass = max(0, biomass)
        substrate = max(0, substrate)
        
        # Bacterial growth and consumption
        mu = self.monod.specific_growth_rate(substrate)
        
        # Biomass change
        dX_dt = mu * biomass
        
        # Substrate uptake by bacteria (mmol/(L·hr))
        substrate_uptake = (mu * biomass) / self.monod.Y_xs
        
        # Bead release rate (mmol/(L·hr))
        bead_release = self.get_bead_release_rate(t)
        
        # Substrate change: dS/dt = bead_release - bacterial_uptake
        dS_dt = bead_release - substrate_uptake
        
        return dX_dt, dS_dt
    
    def simulate(self, initial_od, initial_substrate_mM, time_days):
        """
        Simulate the integrated system.
        
        Parameters:
        -----------
        initial_od : float
            Initial optical density
        initial_substrate_mM : float
            Initial substrate concentration (mM)
        time_days : float
            Simulation time (days)
        
        Returns:
        --------
        dict : Results including time, biomass, substrate, OD
        """
        # Convert OD to biomass
        initial_biomass = self.monod.od_to_biomass(initial_od, self.od_conversion)
        
        # Time points (hourly resolution)
        time_hours = np.linspace(0, time_days * 24, int(time_days * 24 * 4))
        
        # Initial state
        y0 = [initial_biomass, initial_substrate_mM]
        
        # Solve ODE
        solution = odeint(self.integrated_model, y0, time_hours)
        
        biomass = solution[:, 0]
        substrate = solution[:, 1]
        
        # Convert time to days
        time_days_array = time_hours / 24.0
        
        # Convert biomass to OD
        od = self.monod.biomass_to_od(biomass, self.od_conversion)
        
        # Calculate growth rates
        growth_rate = np.array([self.monod.specific_growth_rate(s) for s in substrate])
        
        # Calculate bead release at each time point
        bead_release_rate = np.array([self.get_bead_release_rate(t) for t in time_hours])
        
        # Calculate bacterial uptake at each time point
        bacterial_uptake = np.array([
            (self.monod.specific_growth_rate(s) * b) / self.monod.Y_xs 
            for s, b in zip(substrate, biomass)
        ])
        
        return {
            'time_days': time_days_array,
            'time_hours': time_hours,
            'biomass_g_L': biomass,
            'substrate_mM': substrate,
            'OD600': od,
            'growth_rate_per_hr': growth_rate,
            'bead_release_rate_mmol_L_hr': bead_release_rate,
            'bacterial_uptake_mmol_L_hr': bacterial_uptake,
            'net_substrate_change_mmol_L_hr': bead_release_rate - bacterial_uptake
        }
    
    def plot_results(self, results, save_path=None):
        """
        Plot simulation results.
        
        Parameters:
        -----------
        results : dict
            Results from simulate()
        save_path : str, optional
            Path to save figure
        """
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        time = results['time_days']
        
        # Plot 1: OD over time
        axes[0, 0].plot(time, results['OD600'], 'b-', linewidth=2)
        axes[0, 0].set_xlabel('Time (days)', fontsize=12)
        axes[0, 0].set_ylabel('OD600', fontsize=12)
        axes[0, 0].set_title('Bacterial Growth', fontsize=14, fontweight='bold')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Substrate concentration
        axes[0, 1].plot(time, results['substrate_mM'], 'g-', linewidth=2)
        axes[0, 1].set_xlabel('Time (days)', fontsize=12)
        axes[0, 1].set_ylabel('Formate Concentration (mM)', fontsize=12)
        axes[0, 1].set_title('Substrate Concentration', fontsize=14, fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Plot 3: Growth rate
        axes[1, 0].plot(time, results['growth_rate_per_hr'], 'r-', linewidth=2)
        axes[1, 0].set_xlabel('Time (days)', fontsize=12)
        axes[1, 0].set_ylabel('Specific Growth Rate (hr⁻¹)', fontsize=12)
        axes[1, 0].set_title('Bacterial Growth Rate', fontsize=14, fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Plot 4: Release vs Consumption
        axes[1, 1].plot(time, results['bead_release_rate_mmol_L_hr'], 'b-', 
                       linewidth=2, label='Bead Release')
        axes[1, 1].plot(time, results['bacterial_uptake_mmol_L_hr'], 'r-', 
                       linewidth=2, label='Bacterial Consumption')
        axes[1, 1].plot(time, results['net_substrate_change_mmol_L_hr'], 'g--', 
                       linewidth=2, label='Net Change')
        axes[1, 1].set_xlabel('Time (days)', fontsize=12)
        axes[1, 1].set_ylabel('Rate (mmol/(L·hr))', fontsize=12)
        axes[1, 1].set_title('Substrate Fluxes', fontsize=14, fontweight='bold')
        axes[1, 1].legend(loc='best')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Figure saved to: {save_path}")
        
        plt.show()
    
    def print_summary(self, results):
        """Print summary of simulation results."""
        print("\n" + "="*80)
        print("INTEGRATED SIMULATION SUMMARY")
        print("="*80)
        
        print(f"\nInitial Conditions:")
        print(f"  OD600:              {results['OD600'][0]:.3f}")
        print(f"  Substrate:          {results['substrate_mM'][0]:.2f} mM")
        print(f"  Biomass:            {results['biomass_g_L'][0]:.3f} g/L")
        
        print(f"\nFinal Conditions (Day {results['time_days'][-1]:.1f}):")
        print(f"  OD600:              {results['OD600'][-1]:.3f}")
        print(f"  Substrate:          {results['substrate_mM'][-1]:.2f} mM")
        print(f"  Biomass:            {results['biomass_g_L'][-1]:.3f} g/L")
        print(f"  Growth rate:        {results['growth_rate_per_hr'][-1]:.4f} hr⁻¹")
        
        # Find min/max substrate
        min_substrate = np.min(results['substrate_mM'])
        max_substrate = np.max(results['substrate_mM'])
        min_time = results['time_days'][np.argmin(results['substrate_mM'])]
        max_time = results['time_days'][np.argmax(results['substrate_mM'])]
        
        print(f"\nSubstrate Statistics:")
        print(f"  Minimum:            {min_substrate:.2f} mM (Day {min_time:.1f})")
        print(f"  Maximum:            {max_substrate:.2f} mM (Day {max_time:.1f})")
        print(f"  Average:            {np.mean(results['substrate_mM']):.2f} mM")
        
        # Find max OD
        max_od = np.max(results['OD600'])
        max_od_time = results['time_days'][np.argmax(results['OD600'])]
        
        print(f"\nBacterial Growth:")
        print(f"  Maximum OD:         {max_od:.3f} (Day {max_od_time:.1f})")
        print(f"  Fold increase:      {results['OD600'][-1] / results['OD600'][0]:.1f}x")
        
        if self.beads:
            print(f"\nBead Information:")
            m07_count = sum(1 for b in self.beads if b.bead_type == 'M07')
            m03_count = sum(1 for b in self.beads if b.bead_type == 'M03')
            print(f"  M07 beads:          {m07_count}")
            print(f"  M03 beads:          {m03_count}")
            print(f"  Total beads:        {len(self.beads)}")
        
        print("="*80 + "\n")


def example_integrated_simulation():
    """Example of integrated bead + bacteria simulation."""
    
    print("\n" + "="*80)
    print("INTEGRATED FORMATE BEADS + BACTERIAL CONSUMPTION")
    print("="*80 + "\n")
    
    # Parameters
    volume_ml = 100
    mu_max = 0.3      # 1/hr
    Ks = 2.0          # mM
    Y_xs = 0.02       # g biomass / mmol formate
    
    initial_od = 0.05
    initial_substrate = 1.0  # mM
    simulation_days = 7
    
    # Create integrated model
    model = IntegratedBeadBacteriaModel(volume_ml, mu_max, Ks, Y_xs)
    
    # Calculate bead schedule using ExperimentManager
    print("Step 1: Calculate bead schedule for desired cumulative concentrations...")
    experiment = ExperimentManager(volume_ml)
    
    # Target cumulative concentrations (what we want WITHOUT bacteria)
    desired_cumulative = {
        1: 5.0,
        2: 10.0,
        3: 15.0,
        4: 20.0,
        5: 25.0,
        6: 30.0,
        7: 35.0,
    }
    
    schedule = experiment.calculate_beads_needed(desired_cumulative)
    
    # Set bead schedule in integrated model
    model.set_bead_schedule(schedule)
    
    # Simulate
    print("\nStep 2: Simulate with bacterial consumption...")
    results = model.simulate(initial_od, initial_substrate, simulation_days)
    
    # Print results
    model.print_summary(results)
    
    # Plot results
    print("Generating plots...")
    model.plot_results(results, save_path='integrated_simulation.png')
    
    return model, results


def example_constant_substrate():
    """
    Example: Calculate bead schedule to maintain constant substrate concentration.
    """
    
    print("\n" + "="*80)
    print("EXAMPLE: MAINTAIN CONSTANT SUBSTRATE CONCENTRATION")
    print("="*80 + "\n")
    
    # Parameters
    volume_ml = 100
    mu_max = 0.3          # 1/hr
    Ks = 2.0              # mM
    Y_xs = 0.02           # g/mmol
    target_substrate = 10.0  # mM - target to maintain
    initial_od = 0.05
    days = 7
    
    # Create constant substrate calculator
    calculator = ConstantSubstrateCalculator(
        volume_ml=volume_ml,
        mu_max=mu_max,
        Ks=Ks,
        Y_xs=Y_xs,
        target_substrate_mM=target_substrate
    )
    
    # Estimate OD trajectory (assuming constant substrate)
    od_trajectory = calculator.estimate_od_trajectory(initial_od, days)
    
    print("\nPredicted OD trajectory (with constant substrate):")
    for day, od in od_trajectory.items():
        print(f"  Day {day}: OD = {od:.3f}")
    
    # Calculate bead schedule
    print("\n" + "="*80)
    schedule = calculator.calculate_bead_schedule(
        initial_od=initial_od,
        od_trajectory=od_trajectory,
        days=days
    )
    
    # Now simulate with the calculated schedule
    print("\nStep 2: Verify with simulation...")
    model = IntegratedBeadBacteriaModel(volume_ml, mu_max, Ks, Y_xs)
    model.set_bead_schedule(schedule)
    
    results = model.simulate(initial_od, target_substrate, days)
    
    # Check if substrate stayed near target
    avg_substrate = np.mean(results['substrate_mM'])
    min_substrate = np.min(results['substrate_mM'])
    max_substrate = np.max(results['substrate_mM'])
    
    print("\nVERIFICATION:")
    print(f"  Target substrate:     {target_substrate:.2f} mM")
    print(f"  Average substrate:    {avg_substrate:.2f} mM")
    print(f"  Min substrate:        {min_substrate:.2f} mM")
    print(f"  Max substrate:        {max_substrate:.2f} mM")
    print(f"  Deviation from target: {abs(avg_substrate - target_substrate):.2f} mM")
    
    if abs(avg_substrate - target_substrate) < 2.0:
        print("  ✓ Substrate well-maintained!")
    else:
        print("  ⚠ Consider adjusting bead schedule")
    
    # Plot results
    model.plot_results(results, save_path='constant_substrate_simulation.png')
    
    return calculator, schedule, results


if __name__ == "__main__":
    # Run example with cumulative concentrations
    # model, results = example_integrated_simulation()
    
    # Run example with constant substrate maintenance
    calculator, schedule, results = example_constant_substrate()
