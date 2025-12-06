# Formate Beads + Monod Kinetics Calculator

This project combines formate bead release kinetics with Monod bacterial growth kinetics to predict substrate concentration dynamics in bacterial cultures.

## Overview

The calculator consists of two main components:

1. **Formate Bead Release Calculator** - Calculates how many M07 and M03 beads to add each day to reach target cumulative formate concentrations
2. **Monod Kinetics Calculator** - Models bacterial growth and substrate consumption based on the Monod equation

## Files

### Core Modules

- **`formate_beads_core.py`** - Core logic for bead release calculations
  - Constants: M07 (80mg) and M03 (60mg) beads
  - `Bead` class: Tracks individual bead release and depletion
  - `ExperimentManager` class: Calculates bead schedules for cumulative concentrations

- **`monod_kinetics.py`** - Monod kinetics implementation
  - `MonodKinetics` class: Models bacterial growth and substrate consumption
  - Key equation: μ = μmax × [S] / (Ks + [S])
  - Simulates batch and fed-batch cultures
  - Predicts substrate requirements for target OD

- **`integrated_model.py`** - Combined bead + bacteria model
  - `IntegratedBeadBacteriaModel` class: Simulates both bead release AND bacterial consumption
  - Differential equations for dynamic substrate concentration
  - Generates plots showing growth, substrate, and fluxes

### GUI Applications

- **`formate_beads_gui.py`** - GUI for bead calculations
  - Input: Volume, release profiles, cumulative target concentrations
  - Output: Bead addition schedule (M07 + M03 mix)
  - Shows deviations and bead efficiency

- **`monod_gui.py`** - GUI for Monod kinetics
  - Input: μmax, Ks, Y_x/s, OD, substrate concentration
  - Simulates batch cultures
  - Plots growth curves, substrate consumption
  - Predicts substrate requirements

### Notebooks and Tests

- **`formate_beads_notebook.ipynb`** - Interactive notebook for bead calculations
- **`test_cumulative.py`** - Test script for cumulative concentration tracking
- **`CHANGES.md`** - Documentation of recent updates

## Monod Kinetics Parameters

### Key Parameters

1. **μmax (maximum specific growth rate, hr⁻¹)**
   - Maximum rate at which bacteria can grow
   - Typical values: 0.2-0.6 hr⁻¹ for E. coli
   - Example: 0.3 hr⁻¹ on formate

2. **Ks (half-saturation constant, mM)**
   - Substrate concentration at which μ = μmax/2
   - Lower Ks = higher affinity for substrate
   - Typical values: 0.1-10 mM depending on organism and substrate
   - Example: 2.0 mM for formate

3. **Y_x/s (yield coefficient, g biomass / mmol substrate)**
   - How much biomass is produced per unit substrate consumed
   - Typical values: 0.01-0.05 g/mmol for formate
   - Example: 0.02 g/mmol

4. **OD to Biomass Conversion (g/L per OD unit)**
   - Converts optical density to biomass concentration
   - Typical for E. coli: 0.4 g/L per OD600 unit
   - Can vary by strain and growth conditions

## Usage Examples

### Example 1: Monod Kinetics Only

```python
from monod_kinetics import MonodKinetics

# Create model
monod = MonodKinetics(
    mu_max=0.3,      # hr⁻¹
    Ks=2.0,          # mM
    Y_xs=0.02,       # g/mmol
    volume_ml=100    # mL
)

# Simulate batch culture
results = monod.simulate_batch(
    initial_od=0.1,
    initial_substrate_mM=20.0,
    time_hours=48
)

# Results contain:
# - time_hr: Time points
# - OD600: Optical density
# - substrate_mM: Substrate concentration
# - growth_rate_per_hr: Specific growth rate
# - substrate_consumed_mM: Cumulative consumption
```

### Example 2: Predict Substrate Requirements

```python
from monod_kinetics import MonodKinetics

monod = MonodKinetics(mu_max=0.3, Ks=2.0, Y_xs=0.02, volume_ml=100)

# How much substrate to grow from OD 0.1 to OD 1.0?
prediction = monod.predict_substrate_needed(
    initial_od=0.1,
    final_od=1.0,
    initial_substrate_mM=10.0
)

print(f"Substrate needed: {prediction['substrate_needed_mM']:.2f} mM")
print(f"Sufficient? {prediction['substrate_sufficient']}")
```

### Example 3: Integrated Bead + Bacteria Model

```python
from integrated_model import IntegratedBeadBacteriaModel
from formate_beads_core import ExperimentManager

# Step 1: Calculate bead schedule
experiment = ExperimentManager(volume_ml=100)
schedule = experiment.calculate_beads_needed({
    1: 5.0,   # Day 1: 5 mM cumulative target
    2: 10.0,  # Day 2: 10 mM cumulative target
    3: 15.0,  # ...
    # etc.
})

# Step 2: Create integrated model
model = IntegratedBeadBacteriaModel(
    volume_ml=100,
    mu_max=0.3,
    Ks=2.0,
    Y_xs=0.02
)

# Set bead schedule
model.set_bead_schedule(schedule)

# Step 3: Simulate with bacterial consumption
results = model.simulate(
    initial_od=0.05,
    initial_substrate_mM=1.0,
    time_days=7
)

# Plot results
model.plot_results(results)
```

### Example 4: GUI Applications

```bash
# Run Monod kinetics GUI
python monod_gui.py

# Run formate beads GUI
python formate_beads_gui.py
```

## Understanding the Monod Equation

The Monod equation describes how bacteria grow in response to substrate availability:

```
μ = μmax × [S] / (Ks + [S])
```

Where:
- **μ** = current specific growth rate (hr⁻¹)
- **μmax** = maximum specific growth rate (hr⁻¹)
- **[S]** = substrate concentration (mM)
- **Ks** = half-saturation constant (mM)

### Behavior

1. **High substrate ([S] >> Ks)**: μ ≈ μmax (maximum growth)
2. **Low substrate ([S] << Ks)**: μ ≈ μmax × [S] / Ks (linear with substrate)
3. **At Ks**: μ = μmax / 2 (half-maximum growth)

### Substrate Consumption

Substrate is consumed to produce biomass:

```
dS/dt = -(μ × X) / Y_x/s
```

Where:
- **dS/dt** = rate of substrate consumption
- **X** = biomass concentration (g/L)
- **Y_x/s** = yield coefficient (g biomass / mmol substrate)

### Biomass Growth

Biomass increases according to:

```
dX/dt = μ × X
```

## Integrated Model

The integrated model combines:
1. **Bead release**: Formate is released from beads based on their age and type
2. **Bacterial consumption**: Bacteria consume formate according to Monod kinetics

The net substrate change is:

```
dS/dt = (Bead Release Rate) - (Bacterial Uptake Rate)
```

This allows you to:
- Predict actual substrate concentrations accounting for bacterial consumption
- Optimize bead schedules to maintain target substrate levels
- Understand if beads can keep up with bacterial consumption
- Design fed-batch experiments

## Installation Requirements

```bash
pip install numpy scipy matplotlib tkinter
```

## Tips for Parameter Selection

### μmax (Maximum Growth Rate)
- Literature values for your organism/substrate
- Measure experimentally: ln(OD_final/OD_initial) / time during exponential phase
- Typical E. coli: 0.4-0.6 hr⁻¹ (glucose), 0.2-0.4 hr⁻¹ (other substrates)

### Ks (Half-Saturation Constant)
- Literature values recommended
- Can estimate from growth curves at different substrate concentrations
- Lower Ks = better substrate affinity

### Y_x/s (Yield Coefficient)
- Measure experimentally: ΔBiomass / ΔSubstrate
- Grow culture, measure OD change and substrate consumed
- Convert OD to biomass using conversion factor

### OD Conversion Factor
- Measure experimentally: Dry weight vs OD
- Typical E. coli: 0.3-0.5 g/L per OD600 unit
- Can vary significantly by strain

## Example Workflow

1. **Determine bacterial parameters**
   - Measure μmax, Ks, Y_x/s for your organism
   - Use literature values as starting point

2. **Calculate desired substrate profile**
   - What concentration do you want each day?
   - Use formate beads calculator to get bead schedule

3. **Simulate with bacterial consumption**
   - Use integrated model to see actual substrate dynamics
   - Check if beads can maintain target levels

4. **Optimize**
   - Adjust bead schedule if needed
   - Consider initial substrate concentration
   - Plan for substrate depletion

## Troubleshooting

**Substrate depletes too quickly:**
- Increase bead additions
- Increase initial substrate
- Use more M07 beads (faster release)

**Substrate accumulates (bacteria not consuming):**
- Check if bacteria are growing (low μmax? substrate limitation?)
- Verify Y_x/s is correct (might be overestimating consumption)

**Negative substrate concentrations in simulation:**
- Substrate became limiting
- Increase bead release or initial substrate
- Model will prevent negative values but indicates depletion

## References

- Monod, J. (1949). The growth of bacterial cultures. Annual Review of Microbiology, 3, 371-394.
- Bailey & Ollis (1986). Biochemical Engineering Fundamentals. McGraw-Hill.
- Shuler & Kargi (2002). Bioprocess Engineering: Basic Concepts. Prentice Hall.

## Contact

For questions about this calculator, refer to the code documentation or comments within each module.
