# Changes Summary - Formate Beads Calculator

## Updates Made

### 1. Renamed Bead Types
- **Fast beads** → **M07 beads** (80 mg total formate)
- **Slow beads** → **M03 beads** (60 mg total formate)

### 2. Changed to Cumulative Concentration Mode
**Previous behavior:**
- User specified desired concentration for each day
- Program calculated beads needed to reach that daily concentration
- Formate was assumed to be consumed/cleared each day

**New behavior:**
- User specifies desired **cumulative** concentration for each day
- Formate **accumulates** in the medium over time
- Program tracks total formate and calculates beads needed to reach cumulative target
- Displays both:
  - Today's release (from currently active beads)
  - Current cumulative total (sum of all previous days)

### 3. Files Updated

#### `formate_beads_core.py`
- Renamed constants: `TOTAL_FORMATE_M07_BEAD`, `TOTAL_FORMATE_M03_BEAD`, `TOTAL_MMOL_M07_BEAD`, `TOTAL_MMOL_M03_BEAD`
- Renamed release profiles: `M07_BEAD_RELEASE`, `M03_BEAD_RELEASE`
- Updated `Bead` class to accept 'M07' or 'M03' as bead_type
- Added `cumulative_formate_mmol` tracking to `ExperimentManager`
- Changed method signature: `calculate_beads_needed(desired_cumulative_concentration_mM)`
- Updated calculation loop to:
  - Track cumulative formate in medium
  - Display today's release and cumulative total
  - Calculate beads needed based on cumulative targets
- Updated summary to show final cumulative formate and concentration

#### `formate_beads_gui.py`
- Updated imports to use `M07_BEAD_RELEASE`, `M03_BEAD_RELEASE`
- Changed label: "Strategy: Mixed M07 + M03 Beads (Cumulative Mode)"
- Changed frame title: "Desired Cumulative Concentrations (mM)"
- Updated help text: "Enter desired CUMULATIVE formate concentration for each day:"
- Changed default values to show cumulative progression (5, 10, 15, 20, 25, 30, 35 mM)
- Renamed GUI elements:
  - `fast_entries` → `m07_entries`
  - `slow_entries` → `m03_entries`
  - "Fast Beads" label → "M07 Beads"
  - "Slow Beads" label → "M03 Beads"
- Updated `get_release_profile()` method to use 'M07'/'M03'
- Updated profile updates in `calculate()` method

#### `formate_beads_notebook.ipynb`
- Updated all markdown cells to reference M07/M03 instead of fast/slow
- Updated constants section with M07/M03 naming
- Changed release profile dictionaries to `M07_BEAD_RELEASE`, `M03_BEAD_RELEASE`
- Updated `Bead` class docstring and logic for 'M07'/'M03' types
- Updated `ExperimentManager` class:
  - Added `cumulative_formate_mmol` tracking
  - Changed to cumulative concentration calculation
  - Updated docstrings to explain cumulative mode
- Changed example usage to show cumulative concentrations (5, 10, 15, 20, 25, 30, 35 mM)
- Updated schedule display to show M07/M03 counts
- Updated custom calculation example with cumulative values

## Example Usage

### Input (Cumulative Concentrations):
```python
desired_cumulative_concentration = {
    1: 5.0,   # Day 1: 5 mM total in medium
    2: 10.0,  # Day 2: 10 mM total in medium
    3: 15.0,  # Day 3: 15 mM total in medium
    4: 20.0,  # Day 4: 20 mM total in medium
    5: 25.0,  # Day 5: 25 mM total in medium
    6: 30.0,  # Day 6: 30 mM total in medium
    7: 35.0,  # Day 7: 35 mM total in medium
}
```

### Output Format:
```
Day 1:
  Target cumulative:      5.00 mM (0.500 mmol)
  Today's release:        0.00 mM (0.000 mmol)
  Current cumulative:     0.00 mM (0.000 mmol)
  Active beads:           0 beads
  Additional needed:      5.00 mM (0.500 mmol)
  → ADD X M07 beads + Y M03 beads
  Actual cumulative:      5.XXX mM (0.XXX mmol)
  Deviation:              +/-X.XXX mM (+/-X.XXX mmol)
```

## Key Features Preserved
- ✅ Only whole beads (no fractional beads)
- ✅ Shows deviations from target concentrations
- ✅ Mixed M07/M03 bead strategy
- ✅ Tracks bead depletion and efficiency
- ✅ Editable release profiles
- ✅ GUI, core module, and notebook versions

## Testing Recommendation
Run the program with an example cumulative profile to verify the cumulative tracking works correctly. The cumulative concentration should increase each day as beads continue to release formate into the medium.
