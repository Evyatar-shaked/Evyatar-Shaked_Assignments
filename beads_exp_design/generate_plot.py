"""
Quick plot generator - saves without showing
"""
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

from integrated_model import IntegratedBeadBacteriaModel
import matplotlib.pyplot as plt

print("Generating plot (no display window)...")

# Simple parameters
volume_ml = 100
mu_max = 0.3
Ks = 2.0
Y_xs = 0.02

# Create model
model = IntegratedBeadBacteriaModel(volume_ml, mu_max, Ks, Y_xs)

# Simple bead schedule
schedule = {day: {'M07': 0, 'M03': 1} for day in range(1, 4)}
model.set_bead_schedule(schedule)

# Run simulation
results = model.simulate(initial_od=0.05, initial_substrate_mM=10.0, time_days=3)

# Generate plot (will save automatically)
model.plot_results(results, save_path='quick_plot.png')

# Close the plot to prevent display
plt.close('all')

print(f"\n✓✓✓ SUCCESS! ✓✓✓")
print(f"Plot saved to: quick_plot.png")
print(f"Full path: c:\\Users\\evya1\\OneDrive\\Desktop\\MSc\\Courses\\python\\Evyatar-Shaked_Assignments\\beads_exp_design\\quick_plot.png")
