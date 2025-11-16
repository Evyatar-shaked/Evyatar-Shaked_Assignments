# Lab Dilution Calculator

A professional dilution calculator designed for lab workers to perform C1V1 = C2V2 calculations with ease. Features automatic unit conversion and an intuitive graphical interface.

## 🧪 Features

- **C1V1 = C2V2 Calculations**: Accurate dilution calculations based on the fundamental dilution equation
- **Unit Conversion**: Automatic conversion between different volume and concentration units
- **Multiple Unit Support**:
  - **Volume Units**: L, mL, µL (or uL), nL
  - **Concentration Units**: M, mM, µM (or uM), nM
- **User-Friendly GUI**: Simple interface with dropdown unit selectors
- **Comprehensive Testing**: Full test suite covering business logic and edge cases
- **Clear Instructions**: Step-by-step dilution instructions in results

## 📋 Requirements

- Python 3.6 or higher
- tkinter (usually included with Python)

## 🚀 Installation

### Option 1: Basic Installation (No additional packages needed)

Since tkinter comes with most Python installations, you can run the calculator directly:

```powershell
# Navigate to the day03 directory
cd "c:\Users\evya1\OneDrive\Desktop\MSc\Courses\python\Evyatar-Shaked_Assignments\day03"

# Run the GUI version
python gui_version.py
```

### Option 2: Install Dependencies (if tkinter is missing)

If you get an error about tkinter not being installed:

**On Windows:**
- Tkinter comes pre-installed with Python from python.org
- If missing, reinstall Python with the "tcl/tk and IDLE" option checked

**On Ubuntu/Debian:**
```bash
sudo apt-get install python3-tk
```

**On macOS:**
```bash
brew install python-tk
```

## 📖 Usage

### Using the GUI (Recommended for Lab Work)

1. **Start the application:**
   ```powershell
   python gui_version.py
   ```

2. **Enter your values:**
   - **C1**: Stock concentration (e.g., 10)
   - **C2**: Final desired concentration (e.g., 1)
   - **V2**: Final volume needed (e.g., 100)

3. **Select units** from the dropdown menus for each value

4. **Choose output unit** for the result volumes

5. **Click Calculate** to get your dilution protocol

### Example Calculations

#### Example 1: Simple 10X Dilution
- **Input**: C1 = 10 M, C2 = 1 M, V2 = 100 mL
- **Output**: Use 10 mL stock + 90 mL solvent

#### Example 2: PCR Master Mix
- **Input**: C1 = 10 M, C2 = 1 M, V2 = 50 µL
- **Output**: Use 5 µL master mix + 45 µL water

#### Example 3: Antibody Dilution (1:1000)
- **Input**: C1 = 1000 M, C2 = 1 M, V2 = 10 mL
- **Output**: Use 10 µL antibody + 9990 µL buffer

#### Example 4: Protein Dilution (Mass Units)
- **Input**: C1 = 50 mg/mL, C2 = 5 mg/mL, V2 = 100 mL
- **Output**: Use 10 mL stock + 90 mL buffer

#### Example 5: Drug Dilution with Unit Conversion
- **Input**: C1 = 1 g/mL, C2 = 100 mg/mL, V2 = 500 µL
- **Output**: Use 50 µL stock + 450 µL solvent

### Using as a Python Module

You can also import the functions in your own scripts:

```python
from dilution_core import calculate_dilution, convert_volume, convert_concentration

# Calculate dilution
v1, added = calculate_dilution(
    c1=10,           # Stock concentration
    c2=1,            # Final concentration
    v2=100,          # Final volume
    c1_unit='M',     # Stock concentration unit
    c2_unit='mM',    # Final concentration unit
    v2_unit='mL',    # Final volume unit
    output_unit='µL' # Desired output unit
)

print(f"Use {v1:.2f} µL of stock")
print(f"Add {added:.2f} µL of solvent")

# Convert between units
volume_in_ul = convert_volume(1, 'mL', 'µL')  # 1 mL = 1000 µL
conc_in_um = convert_concentration(1, 'mM', 'µM')  # 1 mM = 1000 µM
```

## 🧪 Running Tests

The project includes comprehensive tests to ensure accuracy of all calculations.

### Run all tests:

```powershell
# Using unittest (no additional packages needed)
python test_dilution.py
```

### Run with pytest (if installed):

```powershell
# Install pytest first (optional)
pip install pytest

# Run tests with pytest
python -m pytest test_dilution.py -v
```

### Test Coverage

The test suite covers:
- ✅ Volume unit conversions (L, mL, µL, nL)
- ✅ Molar concentration conversions (M, mM, µM, nM)
- ✅ Mass concentration conversions (g/mL, mg/mL, µg/mL, ng/mL, g/L, mg/L, µg/L)
- ✅ Basic dilution calculations
- ✅ Mixed unit calculations (within same system)
- ✅ Real-world lab scenarios
- ✅ Edge cases and error handling
- ✅ Validation functions
- ✅ Error handling for mixing molar and mass units

## 📁 Project Structure

```
day03/
├── dilution_core.py      # Core calculation engine with unit conversion
├── gui_version.py        # Graphical user interface
├── test_dilution.py      # Comprehensive test suite
└── README.md            # This file
```

## 🔬 Supported Units

### Volume Units
| Unit | Symbol | Conversion to Liters |
|------|--------|---------------------|
| Liter | L | 1 |
| Milliliter | mL | 0.001 |
| Microliter | µL or uL | 0.000001 |
| Nanoliter | nL | 0.000000001 |

### Concentration Units

#### Molar Concentration Units
| Unit | Symbol | Conversion to Molar |
|------|--------|---------------------|
| Molar | M | 1 |
| Millimolar | mM | 0.001 |
| Micromolar | µM or uM | 0.000001 |
| Nanomolar | nM | 0.000000001 |

#### Mass Concentration Units
| Unit | Symbol | Conversion to g/mL |
|------|--------|-------------------|
| Gram per milliliter | g/mL | 1 |
| Milligram per milliliter | mg/mL | 0.001 |
| Microgram per milliliter | µg/mL or ug/mL | 0.000001 |
| Nanogram per milliliter | ng/mL | 0.000000001 |
| Gram per liter | g/L | 0.001 |
| Milligram per liter | mg/L | 0.000001 |
| Microgram per liter | µg/L or ug/L | 0.000000001 |

**Important:** You can convert between units within the same system (molar to molar, or mass to mass), but **cannot mix** molar and mass concentration units in the same calculation. Both C1 and C2 must use the same unit system.

*Note: Alternative spellings (uL, uM, ug) are supported for systems that don't handle µ character well.*

## 🛠️ API Reference

### `calculate_dilution(c1, c2, v2, c1_unit='M', c2_unit='M', v2_unit='mL', output_unit='mL')`

Calculate dilution volumes based on C1V1 = C2V2.

**Parameters:**
- `c1` (float): Stock concentration
- `c2` (float): Final concentration
- `v2` (float): Final volume
- `c1_unit` (str): Unit for c1 (default: 'M')
- `c2_unit` (str): Unit for c2 (default: 'M')
- `v2_unit` (str): Unit for v2 (default: 'mL')
- `output_unit` (str): Desired unit for output volumes (default: 'mL')

**Returns:**
- `tuple`: (v1, added_volume) both in output_unit

**Raises:**
- `ValueError`: If c1 is zero or units are invalid

### `convert_volume(value, from_unit, to_unit)`

Convert volume from one unit to another.

### `convert_concentration(value, from_unit, to_unit)`

Convert concentration from one unit to another.

### `validate_dilution(c1, c2, v1, v2, tolerance=0.01)`

Validate that a dilution follows C1V1 = C2V2.

## ⚠️ Important Notes

1. **Stock Concentration Must Be Greater Than Final**: For proper dilution, C1 should be greater than C2. If C2 > C1, you're trying to concentrate, not dilute, which will result in negative solvent volumes.

2. **Zero Stock Concentration**: The calculator will raise an error if C1 = 0, as this would result in division by zero.

3. **Unit System Consistency**: Both C1 and C2 must use the same unit system:
   - ✅ **OK**: M to mM, µM to nM, mg/mL to µg/mL
   - ❌ **NOT OK**: M to mg/mL, mM to µg/mL
   - The calculator will raise an error if you try to mix molar and mass units

4. **Unit Conversion**: The calculator automatically handles conversions between units within the same system, so you can use different units for C1 and C2 as long as they're both molar OR both mass-based.

5. **Precision**: All calculations maintain high precision (6+ decimal places) to ensure accuracy for sensitive lab work.

## 🐛 Troubleshooting

### Issue: "tkinter not found" error
**Solution**: Reinstall Python with tkinter support or install python3-tk package (Linux).

### Issue: Calculator gives unexpected results
**Solution**: 
1. Check that C1 > C2 (you're diluting, not concentrating)
2. Verify unit selections are correct
3. Ensure input values are numeric

### Issue: Tests failing
**Solution**: Make sure you're in the correct directory and Python can find all modules:
```powershell
cd "c:\Users\evya1\OneDrive\Desktop\MSc\Courses\python\Evyatar-Shaked_Assignments\day03"
python test_dilution.py
```

### created with claude sonnet 4.5 in the vscode agent mode
### prompt: 
I want to update the program and add more features. The program is designed to help lab workers perform dilutions easily, so I want to include an option for unit conversion (e.g., mL, µL) for both volume and concentration units.

Another improvement is to add a separate test file that verifies the program’s business logic. Additionally, I plan to create a README file that explains the program, including instructions for installing any dependencies if they exist.
