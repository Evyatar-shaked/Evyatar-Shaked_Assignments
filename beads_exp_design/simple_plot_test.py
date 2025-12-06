"""
Simple test to generate a plot file
"""
from integrated_model import IntegratedBeadBacteriaModel
import numpy as np

print("Creating simple simulation...")

# Simple parameters
volume_ml = 100
mu_max = 0.3
Ks = 2.0
Y_xs = 0.02

# Create model
model = IntegratedBeadBacteriaModel(volume_ml, mu_max, Ks, Y_xs)

# Simple bead schedule (1 M03 bead per day)
schedule = {day: {'M07': 0, 'M03': 1} for day in range(1, 4)}
model.set_bead_schedule(schedule)

# Run simulation
print("Running simulation...")
results = model.simulate(initial_od=0.05, initial_substrate_mM=10.0, days=3)

# Generate plot
print("Generating plot...")
model.plot_results(results, save_path='test_plot.png')

print(f"\n✓ Plot saved to: test_plot.png")
print(f"   Location: c:\\Users\\evya1\\OneDrive\\Desktop\\MSc\\Courses\\python\\Evyatar-Shaked_Assignments\\beads_exp_design\\test_plot.png")
