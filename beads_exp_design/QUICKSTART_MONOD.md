# Quick Start Guide - Monod Kinetics Calculator

## Installation

1. Install required packages:
```bash
pip install numpy scipy matplotlib
```

2. Verify tkinter is available (usually comes with Python):
```bash
python -c "import tkinter; print('tkinter OK')"
```

## Quick Examples

### 1. Run Monod Kinetics GUI (Easiest!)

```bash
python monod_gui.py
```

**What it does:**
- Calculate bacterial growth with Monod kinetics
- Predict substrate consumption
- Plot growth curves
- Determine how much substrate you need for target OD

**Input parameters:**
- **μmax**: Maximum growth rate (try 0.3 hr⁻¹)
- **Ks**: Half-saturation constant (try 2.0 mM)
- **Y_x/s**: Yield coefficient (try 0.02 g/mmol)
- **Volume**: Culture volume (try 100 mL)
- **Initial OD**: Starting bacteria density (try 0.1)
- **Initial substrate**: Starting formate (try 20 mM)
- **Time**: How long to simulate (try 48 hours)

### 2. Run Example Script

```bash
python monod_kinetics.py
```

This runs a demo showing:
- Monod parameters summary
- Substrate requirement prediction
- Batch culture simulation

### 3. Run Integrated Model (Advanced)

```bash
python integrated_model.py
```

This simulates both:
- Formate release from beads
- Bacterial consumption of formate

Shows you the actual substrate dynamics when bacteria are present!

### 4. Use in Your Code

```python
from monod_kinetics import MonodKinetics

# Create model
monod = MonodKinetics(
    mu_max=0.3,      # Maximum growth rate (hr⁻¹)
    Ks=2.0,          # Half-saturation constant (mM)
    Y_xs=0.02,       # Yield (g biomass / mmol substrate)
    volume_ml=100    # Culture volume (mL)
)

# Simulate batch culture
results = monod.simulate_batch(
    initial_od=0.1,           # Starting OD600
    initial_substrate_mM=20.0, # Starting formate (mM)
    time_hours=48             # Simulation time
)

# Access results
print(f"Final OD: {results['OD600'][-1]:.2f}")
print(f"Final substrate: {results['substrate_mM'][-1]:.2f} mM")
print(f"Substrate consumed: {results['substrate_consumed_mM'][-1]:.2f} mM")
```

## Understanding Your Results

### Growth Curve (OD vs Time)
- **Lag phase**: Initial flat region (bacteria adapting)
- **Exponential phase**: Steep growth (μ ≈ μmax)
- **Stationary phase**: Plateau (substrate depleted or other limitation)

### Substrate Concentration
- Should decrease over time as bacteria consume it
- When substrate → 0, growth rate → 0 (stationary phase)

### Growth Rate (μ)
- Starts high when substrate is abundant
- Decreases as substrate is consumed
- Follows Monod equation: μ = μmax × [S] / (Ks + [S])

## Common Use Cases

### Case 1: "How much formate do I need?"

Use the **"Predict Substrate Needed"** button in GUI:
- Enter: Initial OD = 0.1, Target OD = 1.0
- Result tells you: Need X mM of formate

### Case 2: "Will my bacteria consume all the formate?"

Run a simulation:
- Set initial substrate to your planned amount
- Look at final substrate concentration
- If > 0, some remains; if ≈ 0, all consumed

### Case 3: "Should I add more beads?"

Use integrated model:
- Set up bead schedule
- Simulate with bacteria
- Check if substrate stays above zero
- If it drops too low, add more beads!

## Typical Parameter Values

### E. coli on Different Substrates

| Substrate | μmax (hr⁻¹) | Ks (mM) | Notes |
|-----------|-------------|---------|-------|
| Glucose   | 0.5-0.7     | 0.1-1.0 | Fast growth |
| Formate   | 0.2-0.4     | 1.0-5.0 | Slower growth |
| Acetate   | 0.3-0.5     | 1.0-3.0 | Medium growth |

### Yield Coefficients (Y_x/s)
- Typical range: 0.01-0.05 g biomass / mmol substrate
- Higher for energy-rich substrates (glucose)
- Lower for simple substrates (formate, acetate)

### OD Conversion
- E. coli: ~0.4 g/L per OD600 unit
- Can vary: 0.3-0.5 typical range
- Measure for your specific strain!

## Troubleshooting

**"Bacteria don't grow in simulation"**
- Check μmax > 0
- Check initial substrate > 0
- Check initial OD > 0

**"Substrate goes negative"**
- This shouldn't happen (model prevents it)
- But indicates severe substrate limitation
- Add more initial substrate or increase bead release

**"Growth is too slow"**
- Increase μmax (check literature values)
- Ensure substrate concentration is above Ks

**"Bacteria grow too much"**
- Check Y_x/s (might be too high)
- Verify substrate consumption makes sense

## Next Steps

1. **Measure your parameters experimentally**
   - Grow culture, measure OD over time → calculate μmax
   - Measure substrate before/after growth → calculate Y_x/s

2. **Use literature values as starting point**
   - Search for your organism + substrate
   - Example: "E. coli formate Monod parameters"

3. **Integrate with bead calculator**
   - Design bead schedule for target concentrations
   - Simulate with bacterial consumption
   - Optimize schedule based on results

## Example Session

```python
# 1. Create model with your parameters
from monod_kinetics import MonodKinetics

monod = MonodKinetics(mu_max=0.3, Ks=2.0, Y_xs=0.02, volume_ml=100)

# 2. Check parameters
monod.print_summary()

# 3. Predict what you need
prediction = monod.predict_substrate_needed(
    initial_od=0.05,
    final_od=1.0,
    initial_substrate_mM=0
)
print(f"Need {prediction['substrate_needed_mM']:.1f} mM formate")

# 4. Simulate to verify
results = monod.simulate_batch(
    initial_od=0.05,
    initial_substrate_mM=prediction['substrate_needed_mM'],
    time_hours=48
)
print(f"Final OD: {results['OD600'][-1]:.2f}")
print(f"Remaining substrate: {results['substrate_mM'][-1]:.2f} mM")
```

## Questions?

See **README_MONOD.md** for detailed documentation!
