# Constant Substrate Concentration Calculator

## What It Does

This calculator determines **how many beads to add each day** to maintain a constant substrate concentration, accounting for:
- **Bead release** - Formate released from M07 and M03 beads
- **Bacterial consumption** - Substrate consumed according to Monod kinetics

The program balances these two rates to keep substrate concentration stable over time.

## Quick Start

### Method 1: Interactive Calculator (Easiest!)

```bash
python constant_substrate_calculator.py
```

You'll be prompted for:
- **μmax**: Maximum growth rate (hr⁻¹) - e.g., 0.3
- **Ks**: Half-saturation constant (mM) - e.g., 2.0
- **Y_x/s**: Yield coefficient (g/mmol) - e.g., 0.02
- **Volume**: Culture volume (mL) - e.g., 100
- **Target substrate**: Concentration to maintain (mM) - e.g., 10.0
- **Initial OD**: Starting bacterial density - e.g., 0.05
- **Days**: Duration of experiment - e.g., 7

The program will:
1. Calculate bead schedule
2. Run simulation to verify
3. Show you a table for your lab notebook

### Method 2: Quick Example

```bash
python constant_substrate_calculator.py --example
```

Runs with default parameters (maintain 10 mM for 7 days).

### Method 3: Python Code

```python
from integrated_model import ConstantSubstrateCalculator

# Create calculator
calculator = ConstantSubstrateCalculator(
    volume_ml=100,
    mu_max=0.3,           # 1/hr
    Ks=2.0,               # mM
    Y_xs=0.02,            # g/mmol
    target_substrate_mM=10.0  # Target concentration
)

# Estimate bacterial growth trajectory
od_trajectory = calculator.estimate_od_trajectory(
    initial_od=0.05,
    time_days=7
)

# Calculate bead schedule
schedule = calculator.calculate_bead_schedule(
    initial_od=0.05,
    od_trajectory=od_trajectory,
    days=7
)

# Schedule is: {1: {'M07': 2, 'M03': 1}, 2: {...}, ...}
```

## How It Works

### The Balance Equation

For each day, the program calculates:

```
Substrate Balance = Bead Release - Bacterial Consumption
```

To maintain constant substrate:
```
Bead Release ≈ Bacterial Consumption
```

### Step-by-Step Process

**Day 1:**
1. Calculate bacterial consumption rate at target concentration
   - Uses Monod equation: μ = μmax × [S] / (Ks + [S])
   - Consumption rate = (μ × biomass) / Y_x/s
2. Check release from existing beads (none on Day 1)
3. Calculate beads needed to match consumption
4. Add M07 and/or M03 beads

**Day 2+:**
1. Calculate consumption (bacteria have grown, so higher consumption)
2. Check release from **all active beads** (including those added on Day 1)
3. Calculate **additional** beads needed
4. Add new beads as needed

### Key Features

✅ **Accounts for bacterial growth** - Consumption increases as OD increases
✅ **Considers bead aging** - Older beads release less (based on release profiles)
✅ **Uses both bead types** - M07 (fast release) and M03 (sustained release)
✅ **Whole beads only** - No fractional beads
✅ **Verification** - Simulates to confirm substrate stays stable

## Understanding the Output

### Example Output

```
Day 1:
  Expected OD:              0.05
  Bacterial consumption:    0.156 mmol/day
  Existing bead release:    0.000 mmol/day
  Additional needed:        0.156 mmol/day
  → ADD 0 M07 beads + 1 M03 beads
  Total release:            0.150 mmol/day
  Net substrate balance:    -0.006 mmol/day
  ⚠ Slight depletion (-0.006 mmol/day)

Day 2:
  Expected OD:              0.10
  Bacterial consumption:    0.312 mmol/day
  Existing bead release:    0.140 mmol/day  ← Bead from Day 1
  Additional needed:        0.172 mmol/day
  → ADD 0 M07 beads + 1 M03 beads
  ...
```

### What to Look For

- **Net substrate balance** should be close to 0
  - **Positive** = substrate accumulating (add fewer beads)
  - **Negative** = substrate depleting (add more beads)
  - **~0** = balanced (perfect!)

- **Total release vs consumption** should match closely

- **Verification results** after simulation:
  - Average substrate should be near target
  - Deviation < 2 mM is good

## Example Scenarios

### Scenario 1: Maintain 10 mM for Exponential Growth

```python
calculator = ConstantSubstrateCalculator(
    volume_ml=100,
    mu_max=0.3,
    Ks=2.0,
    Y_xs=0.02,
    target_substrate_mM=10.0  # Well above Ks
)
```

Result: Bacteria grow at near-maximum rate, high consumption

### Scenario 2: Maintain 1 mM (Limited Growth)

```python
calculator = ConstantSubstrateCalculator(
    volume_ml=100,
    mu_max=0.3,
    Ks=2.0,
    Y_xs=0.02,
    target_substrate_mM=1.0  # Below Ks
)
```

Result: Slower growth, less consumption, fewer beads needed

### Scenario 3: Large Culture Volume

```python
calculator = ConstantSubstrateCalculator(
    volume_ml=500,  # Larger volume
    mu_max=0.3,
    Ks=2.0,
    Y_xs=0.02,
    target_substrate_mM=10.0
)
```

Result: More beads needed (larger total substrate amount)

## Comparison with Other Modes

### Mode 1: Cumulative Concentration (Original)
- User specifies cumulative target (e.g., 5, 10, 15, 20 mM over days)
- Assumes NO bacterial consumption
- Simple bead release calculation

### Mode 2: Constant Concentration (NEW!)
- User specifies constant target (e.g., 10 mM every day)
- Accounts FOR bacterial consumption
- Balances release and consumption
- **More realistic for bacterial cultures**

### When to Use Which?

**Use Cumulative Mode when:**
- No bacteria present
- Testing bead release kinetics
- Abiotic experiments

**Use Constant Mode when:**
- Bacteria are actively growing
- Want steady-state substrate levels
- Designing fed-batch fermentation
- Maintaining optimal growth conditions

## Adjusting the Schedule

If simulation shows substrate is not well-maintained:

### Substrate drops too low?
1. Add more beads per day
2. Use more M07 beads (faster release)
3. Increase target concentration

### Substrate accumulates?
1. Use fewer beads
2. Use more M03 beads (slower, sustained release)
3. Check bacterial parameters (maybe Y_x/s too high?)

### Manual adjustment:
```python
# After getting schedule
schedule[3]['M07'] += 1  # Add one more M07 bead on Day 3
schedule[5]['M03'] -= 1  # Remove one M03 bead from Day 5

# Re-simulate to verify
model.set_bead_schedule(schedule)
results = model.simulate(...)
```

## Tips for Best Results

1. **Measure your parameters accurately**
   - μmax, Ks, Y_x/s are critical
   - Use literature values as starting point
   - Refine with your experimental data

2. **Provide realistic OD trajectory**
   - The calculator can estimate (assuming constant substrate)
   - Or you can provide measured/expected values
   - More accurate OD → better bead schedule

3. **Start conservatively**
   - It's easier to add more beads than remove them
   - Monitor substrate in real experiment
   - Adjust future days based on actual readings

4. **Consider bead preparation time**
   - M07 and M03 beads need to be prepared in advance
   - Calculate schedule before starting experiment
   - Have extra beads ready for adjustments

## Validation

Always simulate after calculating:
```python
# Calculate schedule
schedule = calculator.calculate_bead_schedule(...)

# Verify with simulation
model = IntegratedBeadBacteriaModel(volume_ml, mu_max, Ks, Y_xs)
model.set_bead_schedule(schedule)
results = model.simulate(initial_od, target_substrate, days)

# Check results
avg_substrate = np.mean(results['substrate_mM'])
print(f"Average substrate: {avg_substrate:.2f} mM")
print(f"Target: {target_substrate:.2f} mM")
```

## Troubleshooting

**"Bacteria don't grow in simulation"**
- Check initial OD > 0
- Verify μmax > 0
- Ensure target substrate > 0

**"Substrate drops to zero"**
- Not enough beads
- Consumption rate too high
- Check Y_x/s (might be too low)
- Try higher target concentration

**"Beads accumulate substrate"**
- Too many beads
- Consumption rate too low
- Check μmax and Ks values
- Bacteria might not be growing as expected

**"Schedule requires too many beads"**
- Large culture volume?
- High target concentration?
- Fast-growing bacteria?
- All are legitimate reasons for many beads!

## Advanced Usage

### Custom OD Trajectory

If you have measured/expected OD values:

```python
# Your experimental data or predictions
od_trajectory = {
    1: 0.05,
    2: 0.12,
    3: 0.28,
    4: 0.65,
    5: 1.50,
    6: 2.00,  # Stationary phase
    7: 2.00   # Stationary phase
}

schedule = calculator.calculate_bead_schedule(
    initial_od=0.05,
    od_trajectory=od_trajectory,  # Use your data
    days=7
)
```

### Function-based Trajectory

```python
# Define growth function
def my_od_function(day):
    # Your custom growth model
    if day <= 4:
        return 0.05 * np.exp(0.3 * day)
    else:
        return 2.0  # Stationary phase
    
schedule = calculator.calculate_bead_schedule(
    initial_od=0.05,
    od_trajectory=my_od_function,  # Pass function
    days=7
)
```

## Files Reference

- **constant_substrate_calculator.py** - Interactive calculator (run this!)
- **integrated_model.py** - Contains ConstantSubstrateCalculator class
- **formate_beads_core.py** - Bead release kinetics
- **monod_kinetics.py** - Bacterial growth kinetics

## Next Steps

1. **Run the interactive calculator** with your parameters
2. **Get the bead schedule** for your lab notebook
3. **Verify with simulation** to ensure it works
4. **Prepare your beads** according to the schedule
5. **Monitor substrate** during experiment to validate predictions

## Questions?

See **README_MONOD.md** for more details on Monod kinetics and bacterial parameters.

---

**Summary**: This calculator makes it easy to maintain constant substrate concentration in bacterial cultures by automatically calculating when and how many beads to add, accounting for both bead release kinetics and bacterial consumption dynamics!
