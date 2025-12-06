"""
Visual Guide to Output Plots
"""

print("""
═══════════════════════════════════════════════════════════════════════════════
    OUTPUT PLOTS - VISUAL GUIDE
═══════════════════════════════════════════════════════════════════════════════

When you run the constant substrate calculator with plotting enabled, you get
4 comprehensive plots arranged in a 2x2 grid:


┌─────────────────────────────────────┬─────────────────────────────────────┐
│  PLOT 1: Bacterial Growth           │  PLOT 2: Substrate Concentration    │
│  ─────────────────────────           │  ────────────────────────────        │
│                                      │                                     │
│      OD600 vs Time                   │      Substrate (mM) vs Time         │
│                                      │                                     │
│  2.0 │                    ╱──        │  12 │ ─╲  ╱─╲  ╱─╲                 │
│      │                 ╱─            │     │    ╲╱   ╲╱   ╲                │
│  1.5 │              ╱─               │  10 │─────────────────             │
│      │           ╱─                  │     │                               │
│  1.0 │        ╱─                     │   8 │                               │
│      │     ╱─                        │     │                               │
│  0.5 │  ╱─                           │   6 │                               │
│      │╱                              │     │                               │
│  0.0 └─────────────────────► Days   │   0 └─────────────────────► Days   │
│      0   1   2   3   4   5   6   7  │      0   1   2   3   4   5   6   7 │
│                                      │                                     │
│  Shows: Exponential bacterial        │  Shows: Substrate oscillates        │
│  growth as substrate is maintained   │  around target as beads release     │
│  at constant level                   │  and bacteria consume               │
│                                      │                                     │
└─────────────────────────────────────┴─────────────────────────────────────┘


┌─────────────────────────────────────┬─────────────────────────────────────┐
│  PLOT 3: Growth Rate                 │  PLOT 4: Substrate Fluxes           │
│  ────────────────────                 │  ─────────────────────              │
│                                      │                                     │
│      μ (hr⁻¹) vs Time                │      Rates (mmol/L/hr) vs Time      │
│                                      │                                     │
│ 0.30 │─────────────────              │ 0.04 │                             │
│      │                               │      │  Bead Release ─────          │
│ 0.25 │                               │ 0.03 │  Bacterial Uptake ─────      │
│      │                               │      │  Net Change ─ ─ ─           │
│ 0.20 │                               │ 0.02 │   ╱╲      ╱╲                │
│      │                               │      │  ╱  ╲    ╱  ╲               │
│ 0.15 │                               │ 0.01 │ ╱    ╲  ╱    ╲              │
│      │                               │      │╱      ╲╱      ╲             │
│ 0.10 │╲                              │ 0.00 │────────────────────          │
│      │ ╲                             │      │                             │
│ 0.05 │  ╲                            │-0.01 └─────────────────────► Days  │
│      └─────────────────────► Days   │      0   1   2   3   4   5   6   7 │
│      0   1   2   3   4   5   6   7  │                                     │
│                                      │                                     │
│  Shows: Growth rate stays constant   │  Shows: Balance between bead        │
│  at high substrate, drops as         │  release (input) and bacterial      │
│  substrate depletes                  │  consumption (output)               │
│                                      │                                     │
└─────────────────────────────────────┴─────────────────────────────────────┘


═══════════════════════════════════════════════════════════════════════════════
    WHAT EACH PLOT TELLS YOU
═══════════════════════════════════════════════════════════════════════════════

🔵 PLOT 1 - Bacterial Growth (OD600 vs Time)
   ─────────────────────────────────────────
   • Shows how bacterial population increases over time
   • Should show exponential growth if substrate is well-maintained
   • Slope indicates growth rate (steep = fast growth)
   • Use this to verify bacteria are growing as expected

   What to look for:
   ✓ Smooth exponential curve = good growth
   ✓ Plateaus = growth limitation (substrate too low?)
   ✗ No growth = check parameters or initial conditions


🟢 PLOT 2 - Substrate Concentration (Substrate vs Time)  
   ────────────────────────────────────────────────────
   • Shows substrate concentration over time
   • Should oscillate around your target concentration
   • Dips = bacteria consuming faster than beads release
   • Peaks = beads releasing faster than bacteria consume

   What to look for:
   ✓ Average near target = well-balanced
   ✓ Small oscillations = good control
   ✗ Drops to zero = not enough beads
   ✗ Always increasing = too many beads


🔴 PLOT 3 - Specific Growth Rate (μ vs Time)
   ─────────────────────────────────────────
   • Shows bacterial growth rate (hr⁻¹) over time
   • Directly related to substrate concentration via Monod equation
   • μ = μmax × [S] / (Ks + [S])
   • Should stay relatively constant if substrate is maintained

   What to look for:
   ✓ Stays near μmax = substrate above Ks (good)
   ✓ Constant rate = stable substrate
   ✗ Decreasing = substrate depleting
   ✗ Very low = substrate below Ks


🟡 PLOT 4 - Substrate Fluxes (Rates vs Time)
   ─────────────────────────────────────────
   • Shows three lines:
     - Blue: Bead release rate (input)
     - Red: Bacterial consumption rate (output)
     - Green: Net change (blue - red)
   
   • Net change near zero = balanced (substrate stable)
   • Positive net = substrate accumulating
   • Negative net = substrate depleting

   What to look for:
   ✓ Net near zero = well-balanced schedule
   ✓ Blue and red lines close = good match
   ✗ Large positive net = too many beads
   ✗ Large negative net = not enough beads


═══════════════════════════════════════════════════════════════════════════════
    INTERPRETING RESULTS
═══════════════════════════════════════════════════════════════════════════════

SCENARIO 1: Perfect Balance
────────────────────────────
  Plot 2: Substrate stays at 10 ± 1 mM
  Plot 4: Net flux near zero
  → ✓ Bead schedule is optimal!


SCENARIO 2: Substrate Drops Too Low
────────────────────────────────────
  Plot 2: Substrate decreases over time
  Plot 4: Red line (consumption) > Blue line (release)
  → ⚠ Need more beads! Add more M07 or M03


SCENARIO 3: Substrate Accumulates
──────────────────────────────────
  Plot 2: Substrate increases over time
  Plot 4: Blue line (release) > Red line (consumption)
  → ⚠ Too many beads! Reduce bead additions


SCENARIO 4: Good Average, High Oscillations
────────────────────────────────────────────
  Plot 2: Average is correct but swings ±3 mM
  → ⚠ Consider smoother bead schedule
  → Use more M03 (sustained) vs M07 (burst)


═══════════════════════════════════════════════════════════════════════════════
    ADDITIONAL OUTPUTS
═══════════════════════════════════════════════════════════════════════════════

Besides the plots, you also get:

1. NUMERICAL SUMMARY
   ─────────────────
   • Initial and final OD, substrate, biomass
   • Min/max/average substrate
   • Growth fold-increase
   • Bead counts (M07 and M03)

2. DAILY SCHEDULE
   ──────────────
   Day │ M07 │ M03 │ Notes
   ────┼─────┼─────┼──────────────
    1  │  0  │  1  │ Add to culture
    2  │  0  │  1  │ Add to culture
    3  │  1  │  1  │ Add to culture
   ... 

3. VERIFICATION METRICS
   ────────────────────
   • Target vs Average substrate
   • Deviation from target
   • Success/Warning indicators


═══════════════════════════════════════════════════════════════════════════════
    HOW TO GET THE PLOTS
═══════════════════════════════════════════════════════════════════════════════

METHOD 1: Interactive Calculator
─────────────────────────────────
$ python constant_substrate_calculator.py
...
Run simulation to verify? (y/n): y
...
Generate plots? (y/n): y    ← Say yes here!


METHOD 2: Python Code
──────────────────────
from integrated_model import IntegratedBeadBacteriaModel

model = IntegratedBeadBacteriaModel(...)
model.set_bead_schedule(schedule)
results = model.simulate(...)

# Generate plots
model.plot_results(results, save_path='my_plots.png')


METHOD 3: Demo Script
──────────────────────
$ python demo_constant_substrate.py

(Automatically generates all plots)


═══════════════════════════════════════════════════════════════════════════════
    FILE OUTPUTS
═══════════════════════════════════════════════════════════════════════════════

Plots are saved as high-resolution PNG files:
  • 300 DPI (publication quality)
  • 14" x 10" size
  • All 4 plots in one figure
  • Ready for presentations or papers

Default filenames:
  • constant_substrate_simulation.png
  • demo_constant_substrate.png
  • Or specify your own: model.plot_results(results, save_path='myname.png')


═══════════════════════════════════════════════════════════════════════════════

TIP: Run the demo to see all plots in action!
     $ python demo_constant_substrate.py

═══════════════════════════════════════════════════════════════════════════════
""")
