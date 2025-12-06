"""
Visual Workflow Guide - ASCII Art Edition
"""

print("""
═══════════════════════════════════════════════════════════════════════════════
    FORMATE BEADS + MONOD KINETICS CALCULATOR - WORKFLOW
═══════════════════════════════════════════════════════════════════════════════

SCENARIO 1: Just Calculate Bead Schedule (No Bacteria)
─────────────────────────────────────────────────────────────────────────────
    
    ┌─────────────────────────┐
    │  formate_beads_gui.py   │  ← Run this
    │  or                     │
    │  formate_beads_core.py  │
    └───────────┬─────────────┘
                │
                │ Input: Desired cumulative concentrations
                │        {1: 5mM, 2: 10mM, 3: 15mM...}
                │
                ↓
    ┌───────────────────────────┐
    │  Bead Addition Schedule   │
    │  Day 1: 2 M07 + 1 M03    │
    │  Day 2: 0 M07 + 2 M03    │
    │  Day 3: 1 M07 + 1 M03    │
    │  ...                      │
    └───────────────────────────┘


SCENARIO 2: Understand Bacterial Consumption Only
─────────────────────────────────────────────────────────────────────────────

    ┌─────────────────────────┐
    │   monod_gui.py          │  ← Run this
    │   or                    │
    │   monod_kinetics.py     │
    └───────────┬─────────────┘
                │
                │ Input: μmax, Ks, Y_x/s, initial OD, substrate
                │
                ↓
    ┌───────────────────────────────────────┐
    │  Growth Curves + Substrate Depletion │
    │  • OD vs time                        │
    │  • Substrate vs time                 │
    │  • Growth rate vs time               │
    │  • Total consumption                 │
    └───────────────────────────────────────┘


SCENARIO 3: Integrated Model (Beads + Bacteria) - RECOMMENDED!
─────────────────────────────────────────────────────────────────────────────

    ┌───────────────────────┐
    │ STEP 1: Design Beads  │
    │  formate_beads_core   │
    └──────────┬────────────┘
               │
               │ Calculate schedule
               │
               ↓
    ┌──────────────────────────────┐
    │ Bead Schedule                │
    │ {1: {'M07': 2, 'M03': 1},   │
    │  2: {'M07': 0, 'M03': 2},   │
    │  ...}                        │
    └──────────┬───────────────────┘
               │
               │ Feed schedule to integrated model
               │
               ↓
    ┌───────────────────────────────┐
    │ STEP 2: Run Integrated Model  │
    │  integrated_model.py          │
    │                               │
    │  Combines:                    │
    │  • Bead release kinetics      │
    │  • Bacterial consumption      │
    │  • Monod equation             │
    └──────────┬────────────────────┘
               │
               │ Solve coupled ODEs
               │
               ↓
    ┌─────────────────────────────────────────┐
    │ STEP 3: Analyze Results                │
    │                                         │
    │  Plot 1: OD vs Time                    │
    │  ────────────────────                  │
    │      │     ╱╲                          │
    │  OD  │    ╱  ╲───                      │
    │      │  ╱                              │
    │      └─────────────► Time              │
    │                                         │
    │  Plot 2: Substrate vs Time             │
    │  ────────────────────────              │
    │         │ ╱─╲  ╱─╲                     │
    │  [S]mM  │╱   ╲╱   ╲                    │
    │         │           ╲                  │
    │         └──────────────► Time          │
    │                                         │
    │  Plot 3: Fluxes                        │
    │  ────────────                          │
    │         │ Bead Release ────            │
    │  Rate   │ Bacterial Uptake ────        │
    │         │ Net Change ─ ─ ─            │
    │         └──────────────► Time          │
    │                                         │
    │  Summary:                              │
    │  • Min substrate: X mM (Day Y)         │
    │  • Max OD: Z (Day W)                   │
    │  • Substrate sufficient? Yes/No        │
    └─────────────────────────────────────────┘
               │
               │ If substrate too low:
               │
               ↓
    ┌───────────────────────────┐
    │ STEP 4: Adjust & Iterate  │
    │ • Add more beads          │
    │ • Change bead types       │
    │ • Adjust initial substrate│
    └──────────┬────────────────┘
               │
               │ Go back to STEP 1
               └──────────────────► Repeat until optimal


═══════════════════════════════════════════════════════════════════════════════
    KEY EQUATIONS
═══════════════════════════════════════════════════════════════════════════════

Monod Equation (Bacterial Growth Rate):
    
    μ = μmax × [S] / (Ks + [S])
    
    where:
    • μ = specific growth rate (hr⁻¹)
    • μmax = maximum growth rate (hr⁻¹)
    • [S] = substrate concentration (mM)
    • Ks = half-saturation constant (mM)

Biomass Growth:
    
    dX/dt = μ × X
    
    where:
    • X = biomass concentration (g/L)

Substrate Balance:
    
    dS/dt = R_bead - (μ × X) / Y_x/s
    
    where:
    • R_bead = bead release rate (mmol/(L·hr))
    • Y_x/s = yield coefficient (g biomass / mmol substrate)


═══════════════════════════════════════════════════════════════════════════════
    QUICK REFERENCE - TYPICAL PARAMETER VALUES
═══════════════════════════════════════════════════════════════════════════════

Organism: E. coli
Substrate: Formate

    Parameter               Symbol    Value        Units
    ─────────────────────────────────────────────────────────
    Max growth rate         μmax      0.2-0.4      hr⁻¹
    Half-saturation const   Ks        1.0-5.0      mM
    Yield coefficient       Y_x/s     0.01-0.03    g/mmol
    OD conversion           -         0.4          g/L per OD
    Doubling time           td        1.7-3.5      hours
    
    Note: Measure these experimentally for your specific strain!


═══════════════════════════════════════════════════════════════════════════════
    FILE GUIDE
═══════════════════════════════════════════════════════════════════════════════

Core Modules:
    formate_beads_core.py    - Bead release calculations
    monod_kinetics.py        - Monod kinetics implementation  
    integrated_model.py      - Combined beads + bacteria

GUI Applications:
    formate_beads_gui.py     - Bead calculator GUI
    monod_gui.py            - Monod calculator GUI

Documentation:
    README_MONOD.md          - Complete documentation
    QUICKSTART_MONOD.md      - Quick start guide
    MONOD_SUMMARY.md         - This summary
    CHANGES.md              - Update history

Tests & Examples:
    test_cumulative.py       - Test cumulative tracking
    formate_beads_notebook.ipynb - Interactive notebook


═══════════════════════════════════════════════════════════════════════════════
    GETTING STARTED - 3 EASY STEPS
═══════════════════════════════════════════════════════════════════════════════

Step 1: Install dependencies
    $ pip install numpy scipy matplotlib

Step 2: Run a GUI
    $ python monod_gui.py
    
    OR
    
    $ python formate_beads_gui.py

Step 3: Run integrated example
    $ python integrated_model.py


═══════════════════════════════════════════════════════════════════════════════
    TROUBLESHOOTING
═══════════════════════════════════════════════════════════════════════════════

Problem: "Substrate depletes too quickly"
    Solution: 
    • Increase bead additions
    • Use more M07 beads (faster release)
    • Increase initial substrate

Problem: "Bacteria don't grow"
    Solution:
    • Check μmax > 0
    • Ensure initial substrate > Ks
    • Verify Y_x/s is reasonable

Problem: "Substrate accumulates (not consumed)"
    Solution:
    • Check bacterial parameters (μmax, Ks)
    • Verify initial OD > 0
    • May need lower Y_x/s value


═══════════════════════════════════════════════════════════════════════════════

Made with ♥ for bacterial fermentation experiments!

For questions, see README_MONOD.md or QUICKSTART_MONOD.md

═══════════════════════════════════════════════════════════════════════════════
""")
