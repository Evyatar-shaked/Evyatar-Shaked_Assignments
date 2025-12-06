"""
Demo: Constant Substrate Maintenance with Plots

This script demonstrates the plotting capabilities.
"""

from integrated_model import ConstantSubstrateCalculator, IntegratedBeadBacteriaModel
import numpy as np

print("\n" + "="*80)
print("DEMO: Maintain 10 mM Substrate Over 7 Days")
print("="*80 + "\n")

# Setup (using more realistic parameters for visualization)
volume_ml = 100
mu_max = 0.3
Ks = 2.0
Y_xs = 0.02
target_substrate = 10.0
initial_od = 0.05
days = 3  # Shorter timeframe to avoid extreme growth

# Step 1: Calculate bead schedule
print("STEP 1: Calculate Bead Schedule")
print("-" * 80)
calculator = ConstantSubstrateCalculator(
    volume_ml=volume_ml,
    mu_max=mu_max,
    Ks=Ks,
    Y_xs=Y_xs,
    target_substrate_mM=target_substrate
)

od_trajectory = calculator.estimate_od_trajectory(initial_od, days)
schedule = calculator.calculate_bead_schedule(
    initial_od=initial_od,
    od_trajectory=od_trajectory,
    days=days
)

# Step 2: Run simulation
print("\nSTEP 2: Run Simulation with Plots")
print("-" * 80)
model = IntegratedBeadBacteriaModel(volume_ml, mu_max, Ks, Y_xs)
model.set_bead_schedule(schedule)

results = model.simulate(initial_od, target_substrate, days)

# Step 3: Show statistics
print("\nSTEP 3: Results Summary")
print("-" * 80)
avg_substrate = np.mean(results['substrate_mM'])
min_substrate = np.min(results['substrate_mM'])
max_substrate = np.max(results['substrate_mM'])
final_od = results['OD600'][-1]

print(f"\nTarget Substrate:        {target_substrate:.2f} mM")
print(f"Average Substrate:       {avg_substrate:.2f} mM")
print(f"Min Substrate:           {min_substrate:.2f} mM")
print(f"Max Substrate:           {max_substrate:.2f} mM")
print(f"Deviation:               {abs(avg_substrate - target_substrate):.2f} mM\n")

print(f"Initial OD:              {initial_od:.3f}")
print(f"Final OD:                {final_od:.3f}")
print(f"Growth:                  {final_od / initial_od:.1f}x increase\n")

if abs(avg_substrate - target_substrate) < 2.0:
    print("✓ SUCCESS: Substrate concentration well-maintained!")
else:
    print("⚠ Warning: Substrate deviated from target")

# Step 4: Generate plots
print("\nSTEP 4: Generating Plots...")
print("-" * 80)
print("\nThe following 4 plots will be displayed:")
print("  1. OD600 vs Time (Bacterial Growth Curve)")
print("  2. Substrate Concentration vs Time (Substrate Curve)")
print("  3. Specific Growth Rate vs Time")
print("  4. Substrate Fluxes (Release vs Consumption)\n")

model.plot_results(results, save_path='demo_constant_substrate.png')

print("\n✓ Plots saved to: demo_constant_substrate.png")
print("\n" + "="*80)
print("DEMO COMPLETE!")
print("="*80)
print("\nThe plots show:")
print("  • How bacteria grow over time (OD curve)")
print("  • How substrate concentration is maintained (substrate curve)")
print("  • The balance between bead release and bacterial consumption")
print("\n" + "="*80 + "\n")
