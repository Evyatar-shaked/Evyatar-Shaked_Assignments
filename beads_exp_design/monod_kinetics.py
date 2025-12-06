"""
Monod Kinetics Calculator

This module implements Monod kinetics to model bacterial growth and substrate consumption.
The Monod equation describes the relationship between substrate concentration and 
specific growth rate of bacteria.

Monod equation: μ = μmax * [S] / (Ks + [S])

Where:
- μ = specific growth rate (1/hr)
- μmax = maximum specific growth rate (1/hr)
- [S] = substrate concentration (mM)
- Ks = half-saturation constant (mM) - substrate concentration at μ = μmax/2
"""

import numpy as np
from scipy.integrate import odeint


class MonodKinetics:
    """
    Models bacterial growth and substrate consumption using Monod kinetics.
    """
    
    def __init__(self, mu_max, Ks, Y_xs, volume_ml):
        """
        Initialize Monod kinetics parameters.
        
        Parameters:
        -----------
        mu_max : float
            Maximum specific growth rate (1/hr)
        Ks : float
            Half-saturation constant (mM) - substrate concentration at μ = μmax/2
        Y_xs : float
            Yield coefficient (g biomass / mmol substrate)
            How much biomass is produced per unit substrate consumed
        volume_ml : float
            Culture volume (mL)
        """
        self.mu_max = mu_max
        self.Ks = Ks
        self.Y_xs = Y_xs
        self.volume_ml = volume_ml
        self.volume_L = volume_ml / 1000.0
    
    def od_to_biomass(self, od, conversion_factor=0.4):
        """
        Convert optical density (OD600) to biomass concentration.
        
        Parameters:
        -----------
        od : float
            Optical density at 600 nm
        conversion_factor : float
            Conversion factor (g/L per OD unit)
            Default 0.4 g/L per OD (typical for E. coli)
        
        Returns:
        --------
        float : Biomass concentration (g/L)
        """
        return od * conversion_factor
    
    def biomass_to_od(self, biomass, conversion_factor=0.4):
        """
        Convert biomass concentration to optical density.
        
        Parameters:
        -----------
        biomass : float
            Biomass concentration (g/L)
        conversion_factor : float
            Conversion factor (g/L per OD unit)
        
        Returns:
        --------
        float : Optical density at 600 nm
        """
        return biomass / conversion_factor
    
    def specific_growth_rate(self, substrate_conc):
        """
        Calculate specific growth rate using Monod equation.
        
        Parameters:
        -----------
        substrate_conc : float
            Substrate concentration (mM)
        
        Returns:
        --------
        float : Specific growth rate (1/hr)
        """
        if substrate_conc < 0:
            return 0.0
        return self.mu_max * substrate_conc / (self.Ks + substrate_conc)
    
    def substrate_uptake_rate(self, biomass, substrate_conc):
        """
        Calculate substrate uptake rate.
        
        Parameters:
        -----------
        biomass : float
            Biomass concentration (g/L)
        substrate_conc : float
            Substrate concentration (mM)
        
        Returns:
        --------
        float : Substrate uptake rate (mmol/(L·hr))
        """
        mu = self.specific_growth_rate(substrate_conc)
        # Substrate uptake = (μ * X) / Y_xs
        return (mu * biomass) / self.Y_xs
    
    def monod_model(self, state, t, substrate_addition_rate=0):
        """
        Differential equations for Monod kinetics with continuous substrate addition.
        
        Parameters:
        -----------
        state : tuple
            (biomass, substrate) concentrations
        t : float
            Time (hours)
        substrate_addition_rate : float
            Rate of substrate addition (mmol/(L·hr))
        
        Returns:
        --------
        tuple : (dX/dt, dS/dt)
        """
        biomass, substrate = state
        
        # Ensure non-negative values
        biomass = max(0, biomass)
        substrate = max(0, substrate)
        
        # Calculate specific growth rate
        mu = self.specific_growth_rate(substrate)
        
        # Biomass change: dX/dt = μ * X
        dX_dt = mu * biomass
        
        # Substrate change: dS/dt = -(μ * X) / Y_xs + substrate_addition
        dS_dt = -(mu * biomass) / self.Y_xs + substrate_addition_rate
        
        return dX_dt, dS_dt
    
    def simulate_batch(self, initial_od, initial_substrate_mM, time_hours):
        """
        Simulate batch culture (no substrate addition).
        
        Parameters:
        -----------
        initial_od : float
            Initial optical density
        initial_substrate_mM : float
            Initial substrate concentration (mM)
        time_hours : float or array
            Time points to simulate (hours)
        
        Returns:
        --------
        dict : Contains time, biomass, substrate, OD, and growth_rate arrays
        """
        # Convert OD to biomass
        initial_biomass = self.od_to_biomass(initial_od)
        
        # Create time points
        if isinstance(time_hours, (int, float)):
            t = np.linspace(0, time_hours, 100)
        else:
            t = np.array(time_hours)
        
        # Initial state
        y0 = [initial_biomass, initial_substrate_mM]
        
        # Solve ODE
        solution = odeint(self.monod_model, y0, t)
        
        biomass = solution[:, 0]
        substrate = solution[:, 1]
        
        # Convert biomass to OD
        od = self.biomass_to_od(biomass)
        
        # Calculate growth rates
        growth_rate = np.array([self.specific_growth_rate(s) for s in substrate])
        
        return {
            'time_hr': t,
            'biomass_g_L': biomass,
            'substrate_mM': substrate,
            'OD600': od,
            'growth_rate_per_hr': growth_rate,
            'substrate_consumed_mM': initial_substrate_mM - substrate
        }
    
    def simulate_fed_batch(self, initial_od, initial_substrate_mM, 
                          time_hours, feed_schedule):
        """
        Simulate fed-batch culture with substrate feeding.
        
        Parameters:
        -----------
        initial_od : float
            Initial optical density
        initial_substrate_mM : float
            Initial substrate concentration (mM)
        time_hours : float
            Total simulation time (hours)
        feed_schedule : dict
            Dictionary mapping time (hr) to substrate addition (mM)
            Example: {0: 5.0, 24: 5.0, 48: 5.0}
        
        Returns:
        --------
        dict : Contains time, biomass, substrate, OD, and growth_rate arrays
        """
        # Convert OD to biomass
        initial_biomass = self.od_to_biomass(initial_od)
        
        # Sort feed times
        feed_times = sorted(feed_schedule.keys())
        
        # Results arrays
        all_time = []
        all_biomass = []
        all_substrate = []
        
        current_biomass = initial_biomass
        current_substrate = initial_substrate_mM
        current_time = 0
        
        # Simulate between feeding events
        for i, feed_time in enumerate(feed_times + [time_hours]):
            if feed_time <= current_time:
                continue
            
            # Time points for this interval
            t_interval = np.linspace(current_time, feed_time, 50)
            
            # Solve for this interval
            y0 = [current_biomass, current_substrate]
            solution = odeint(self.monod_model, y0, t_interval)
            
            # Store results
            all_time.extend(t_interval)
            all_biomass.extend(solution[:, 0])
            all_substrate.extend(solution[:, 1])
            
            # Update current state
            current_time = feed_time
            current_biomass = solution[-1, 0]
            current_substrate = solution[-1, 1]
            
            # Add substrate if this is a feed time
            if feed_time in feed_schedule:
                current_substrate += feed_schedule[feed_time]
        
        # Convert to arrays
        time_array = np.array(all_time)
        biomass_array = np.array(all_biomass)
        substrate_array = np.array(all_substrate)
        
        # Convert biomass to OD
        od_array = self.biomass_to_od(biomass_array)
        
        # Calculate growth rates
        growth_rate = np.array([self.specific_growth_rate(s) for s in substrate_array])
        
        return {
            'time_hr': time_array,
            'biomass_g_L': biomass_array,
            'substrate_mM': substrate_array,
            'OD600': od_array,
            'growth_rate_per_hr': growth_rate,
            'feed_times': feed_times,
            'feed_amounts': [feed_schedule[t] for t in feed_times]
        }
    
    def predict_substrate_needed(self, initial_od, final_od, initial_substrate_mM=0):
        """
        Predict how much substrate is needed to grow from initial_od to final_od.
        
        Parameters:
        -----------
        initial_od : float
            Starting optical density
        final_od : float
            Target optical density
        initial_substrate_mM : float
            Initial substrate concentration (mM)
        
        Returns:
        --------
        dict : Contains substrate needed, final substrate, and growth info
        """
        initial_biomass = self.od_to_biomass(initial_od)
        final_biomass = self.od_to_biomass(final_od)
        
        # Biomass change
        delta_biomass = final_biomass - initial_biomass  # g/L
        
        # Substrate needed (mmol/L)
        substrate_needed = delta_biomass / self.Y_xs
        
        # Final substrate concentration
        final_substrate = initial_substrate_mM - substrate_needed
        
        return {
            'initial_OD': initial_od,
            'final_OD': final_od,
            'initial_biomass_g_L': initial_biomass,
            'final_biomass_g_L': final_biomass,
            'biomass_increase_g_L': delta_biomass,
            'substrate_needed_mM': substrate_needed,
            'initial_substrate_mM': initial_substrate_mM,
            'final_substrate_mM': final_substrate,
            'substrate_sufficient': final_substrate >= 0
        }
    
    def print_summary(self):
        """Print a summary of the kinetic parameters."""
        print("\n" + "="*70)
        print("MONOD KINETICS PARAMETERS")
        print("="*70)
        print(f"Maximum specific growth rate (μmax): {self.mu_max:.3f} hr⁻¹")
        print(f"Half-saturation constant (Ks):       {self.Ks:.2f} mM")
        print(f"Yield coefficient (Y_x/s):           {self.Y_xs:.3f} g biomass/mmol substrate")
        print(f"Culture volume:                      {self.volume_ml:.1f} mL ({self.volume_L:.3f} L)")
        print(f"\nAt Ks concentration ({self.Ks} mM):")
        print(f"  Growth rate = {self.mu_max/2:.3f} hr⁻¹ (50% of maximum)")
        print(f"\nDoubling time at μmax:               {np.log(2)/self.mu_max:.2f} hours")
        print("="*70 + "\n")


def example_usage():
    """Example usage of Monod kinetics."""
    
    print("\n" + "="*70)
    print("MONOD KINETICS EXAMPLE - E. coli on Formate")
    print("="*70 + "\n")
    
    # Example parameters for E. coli on formate
    mu_max = 0.3      # 1/hr (typical for E. coli on formate)
    Ks = 2.0          # mM (half-saturation constant)
    Y_xs = 0.02       # g biomass / mmol formate (example value)
    volume = 100      # mL
    
    # Create Monod model
    monod = MonodKinetics(mu_max, Ks, Y_xs, volume)
    monod.print_summary()
    
    # Example 1: Predict substrate needed
    print("EXAMPLE 1: Substrate Requirement")
    print("-" * 70)
    prediction = monod.predict_substrate_needed(
        initial_od=0.1,
        final_od=1.0,
        initial_substrate_mM=10.0
    )
    
    print(f"To grow from OD {prediction['initial_OD']} to OD {prediction['final_OD']}:")
    print(f"  Biomass increase:        {prediction['biomass_increase_g_L']:.3f} g/L")
    print(f"  Substrate needed:        {prediction['substrate_needed_mM']:.2f} mM")
    print(f"  Initial substrate:       {prediction['initial_substrate_mM']:.2f} mM")
    print(f"  Final substrate:         {prediction['final_substrate_mM']:.2f} mM")
    print(f"  Substrate sufficient:    {prediction['substrate_sufficient']}")
    
    # Example 2: Batch simulation
    print("\n\nEXAMPLE 2: Batch Culture Simulation")
    print("-" * 70)
    results = monod.simulate_batch(
        initial_od=0.1,
        initial_substrate_mM=20.0,
        time_hours=48
    )
    
    print(f"Initial conditions:")
    print(f"  OD600:              {results['OD600'][0]:.3f}")
    print(f"  Substrate:          {results['substrate_mM'][0]:.2f} mM")
    
    print(f"\nFinal conditions (after {results['time_hr'][-1]:.1f} hours):")
    print(f"  OD600:              {results['OD600'][-1]:.3f}")
    print(f"  Substrate:          {results['substrate_mM'][-1]:.2f} mM")
    print(f"  Substrate consumed: {results['substrate_consumed_mM'][-1]:.2f} mM")
    print(f"  Final growth rate:  {results['growth_rate_per_hr'][-1]:.4f} hr⁻¹")
    
    # Find maximum OD
    max_od_idx = np.argmax(results['OD600'])
    print(f"\nMaximum OD reached:")
    print(f"  Time:               {results['time_hr'][max_od_idx]:.1f} hours")
    print(f"  OD600:              {results['OD600'][max_od_idx]:.3f}")
    
    print("\n" + "="*70 + "\n")


if __name__ == "__main__":
    example_usage()
