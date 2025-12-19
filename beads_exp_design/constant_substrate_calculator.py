"""
Constant Substrate Calculator

This script calculates the bead schedule needed to maintain a constant
substrate concentration despite bacterial consumption.

USAGE:
------
python constant_substrate_calculator.py

Then enter your parameters when prompted.
"""

from integrated_model import ConstantSubstrateCalculator, IntegratedBeadBacteriaModel
import numpy as np


def interactive_calculator():
    """Interactive calculator for constant substrate maintenance."""
    
    print("\n" + "="*80)
    print("CONSTANT SUBSTRATE CONCENTRATION - BEAD CALCULATOR")
    print("="*80)
    print("\nThis calculator determines how many beads to add each day")
    print("to maintain a constant substrate concentration as bacteria consume it.")
    print("\n" + "="*80 + "\n")
    
    # Get user inputs
    print("BACTERIAL PARAMETERS:")
    print("-" * 40)
    
    mu_max = float(input("μmax (maximum growth rate, hr⁻¹) [default 0.3]: ") or "0.3")
    Ks = float(input("Ks (half-saturation constant, mM) [default 2.0]: ") or "2.0")
    Y_xs = float(input("Y_x/s (yield, g biomass/mmol substrate) [default 0.02]: ") or "0.02")
    
    print("\nEXPERIMENT PARAMETERS:")
    print("-" * 40)
    
    volume_ml = float(input("Culture volume (mL) [default 100]: ") or "100")
    target_substrate = float(input("Target substrate concentration (mM) [default 10.0]: ") or "10.0")
    initial_od = float(input("Initial OD600 [default 0.05]: ") or "0.05")
    days = int(input("Number of days [default 7]: ") or "7")
    
    print("\n" + "="*80)
    print("CALCULATING BEAD SCHEDULE...")
    print("="*80 + "\n")
    
    # Create calculator
    calculator = ConstantSubstrateCalculator(
        volume_ml=volume_ml,
        mu_max=mu_max,
        Ks=Ks,
        Y_xs=Y_xs,
        target_substrate_mM=target_substrate
    )
    
    # Estimate OD trajectory
    od_trajectory = calculator.estimate_od_trajectory(initial_od, days)
    
    # Calculate bead schedule
    schedule = calculator.calculate_bead_schedule(
        initial_od=initial_od,
        od_trajectory=od_trajectory,
        days=days
    )
    
    # Ask if user wants to simulate
    print("\n" + "="*80)
    simulate = input("\nRun simulation to verify? (y/n) [default y]: ") or "y"
    
    if simulate.lower() == 'y':
        print("\nRUNNING SIMULATION...")
        print("="*80 + "\n")
        
        model = IntegratedBeadBacteriaModel(volume_ml, mu_max, Ks, Y_xs)
        model.set_bead_schedule(schedule)
        
        results = model.simulate(initial_od, target_substrate, days)
        
        # Print verification
        avg_substrate = np.mean(results['substrate_mM'])
        min_substrate = np.min(results['substrate_mM'])
        max_substrate = np.max(results['substrate_mM'])
        
        print("\n" + "="*80)
        print("SIMULATION RESULTS")
        print("="*80)
        print(f"\nTarget substrate:        {target_substrate:.2f} mM")
        print(f"Average substrate:       {avg_substrate:.2f} mM")
        print(f"Min substrate:           {min_substrate:.2f} mM")
        print(f"Max substrate:           {max_substrate:.2f} mM")
        print(f"Deviation from target:   {abs(avg_substrate - target_substrate):.2f} mM")
        
        if abs(avg_substrate - target_substrate) < 2.0:
            print("\n✓ Substrate concentration well-maintained!")
        elif avg_substrate < target_substrate - 2.0:
            print("\n⚠ Substrate drops below target - consider adding more beads")
        else:
            print("\n⚠ Substrate accumulates above target - consider fewer beads")
        
        final_od = results['OD600'][-1]
        print(f"\nFinal OD600:             {final_od:.3f}")
        print(f"Fold increase:           {final_od / initial_od:.1f}x")
        
        # Ask about plotting
        plot = input("\nGenerate plots? (y/n) [default y]: ") or "y"
        if plot.lower() == 'y':
            model.plot_results(results, save_path='constant_substrate_results.png')
            print("\n✓ Plots saved to: constant_substrate_results.png")
        
        print("\n" + "="*80)
    
    # Print schedule summary for lab use
    print("\n" + "="*80)
    print("BEAD ADDITION SCHEDULE (FOR LAB NOTEBOOK)")
    print("="*80)
    print(f"\nTarget: Maintain {target_substrate:.2f} mM formate")
    print(f"Volume: {volume_ml} mL")
    print(f"Duration: {days} days\n")
    
    print("Day │ M07 Beads │ M03 Beads │ Notes")
    print("────┼───────────┼───────────┼─────────────────")
    for day in range(1, days + 1):
        m07 = schedule[day]['M07']
        m03 = schedule[day]['M03']
        if m07 > 0 or m03 > 0:
            print(f" {day:2d} │    {m07:2d}     │    {m03:2d}     │ Add to culture")
        else:
            print(f" {day:2d} │     -     │     -     │ No addition needed")
    
    total_m07 = sum(day['M07'] for day in schedule.values())
    total_m03 = sum(day['M03'] for day in schedule.values())
    print("────┴───────────┴───────────┴─────────────────")
    print(f"Total: {total_m07} M07 beads, {total_m03} M03 beads\n")
    
    print("="*80)
    print("\nDone! Use this schedule to maintain constant substrate concentration.")
    print("="*80 + "\n")


def quick_example():
    """Quick example with default parameters."""
    
    print("\n" + "="*80)
    print("QUICK EXAMPLE: Maintain 10 mM substrate for 7 days")
    print("="*80 + "\n")
    
    calculator = ConstantSubstrateCalculator(
        volume_ml=100,
        mu_max=0.3,
        Ks=2.0,
        Y_xs=0.02,
        target_substrate_mM=10.0
    )
    
    od_trajectory = calculator.estimate_od_trajectory(initial_od=0.05, time_days=7)
    
    schedule = calculator.calculate_bead_schedule(
        initial_od=0.05,
        od_trajectory=od_trajectory,
        days=7
    )
    
    return calculator, schedule


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1 and sys.argv[1] == "--example":
        # Run quick example
        quick_example()
    else:
        # Run interactive calculator
        interactive_calculator()
