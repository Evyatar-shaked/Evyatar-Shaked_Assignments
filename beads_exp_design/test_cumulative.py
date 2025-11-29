"""
Test script to verify cumulative concentration tracking
"""

from formate_beads_core import ExperimentManager

# Test with 100 mL volume and cumulative concentrations
volume = 100  # mL

# Cumulative concentration targets (formate accumulates over time)
cumulative_targets = {
    1: 5.0,   # Day 1: 5 mM cumulative
    2: 10.0,  # Day 2: 10 mM cumulative
    3: 15.0,  # Day 3: 15 mM cumulative
    4: 20.0,  # Day 4: 20 mM cumulative
    5: 25.0,  # Day 5: 25 mM cumulative
    6: 30.0,  # Day 6: 30 mM cumulative
    7: 35.0,  # Day 7: 35 mM cumulative
}

print("="*80)
print("TESTING CUMULATIVE CONCENTRATION MODE")
print("="*80)
print(f"\nVolume: {volume} mL")
print("\nCumulative targets:")
for day, conc in cumulative_targets.items():
    print(f"  Day {day}: {conc} mM")

# Create experiment and calculate
experiment = ExperimentManager(volume_ml=volume)
schedule = experiment.calculate_beads_needed(cumulative_targets)

print("\n" + "="*80)
print("VERIFICATION")
print("="*80)
print(f"\nFinal cumulative formate: {experiment.cumulative_formate_mmol:.3f} mmol")
print(f"Final concentration: {experiment.cumulative_formate_mmol / experiment.volume_L:.2f} mM")
print(f"Target was: {cumulative_targets[7]:.2f} mM")
print(f"Difference: {(experiment.cumulative_formate_mmol / experiment.volume_L) - cumulative_targets[7]:.3f} mM")
