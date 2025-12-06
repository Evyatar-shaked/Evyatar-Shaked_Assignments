# Lab Dilution Calculator

A professional dilution calculator designed for lab workers to perform C1V1 = C2V2 calculations with ease. Features automatic unit conversion and an intuitive graphical interface.

## 🧪 Features

- **C1V1 = C2V2 Calculations**: Accurate dilution calculations based on the fundamental dilution equation
- **Unit Conversion**: Automatic conversion between different volume and concentration units
- **Multiple Unit Support**:
  - **Volume Units**: L, mL, µL (or uL), nL
  - **Concentration Units**: M, mM, µM (or uM), nM, mg/mL, µg/mL, g/L
- **Batch Processing**: Process multiple dilutions from Excel files at once
- **User-Friendly GUI**: Simple interface with dropdown unit selectors
- **Comprehensive Testing**: Full test suite covering business logic and edge cases
- **Clear Instructions**: Step-by-step dilution instructions in results

## 📋 Requirements

- Python 3.6 or higher
- tkinter (usually included with Python)
- pandas and openpyxl (for Excel batch processing)

## 🚀 Installation

### Option 1: GUI Only (Basic Installation)

Since tkinter comes with most Python installations, you can run the calculator directly:

```powershell
# Run the GUI version
python gui_version.py
```

### Option 2: Full Installation (GUI + Excel Batch Processing)

For Excel batch processing capabilities, install the required packages:

```powershell
# Install all dependencies
pip install -r requirements.txt

# Or install individually
pip install pandas openpyxl
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

### Batch Processing from Excel Files

For processing multiple dilutions at once, you can use an Excel spreadsheet.

#### Excel File Format (Recommended: Units in Headers)

Create an Excel file (.xlsx) with units specified in the column headers. This is the recommended format as it's cleaner and more typical for lab work:

| C1 (M) | C2 (mM) | V2 (mL) |
|--------|---------|---------|
| 10 | 1000 | 100 |
| 2 | 500 | 200 |
| 0.1 | 10 | 50 |
| 5 | 1000 | 10 |

**Header format options:**
- `C1 (M)` - Parentheses format (recommended)
- `C1_M` - Underscore format
- `C1-M` - Dash format
- `C1 M` - Space format

**Important rules:**
- Specify units in column headers using one of the formats above
- Cells contain only numeric values (e.g., `10`, `100`, `0.5`)
- Base column names must be: `C1`, `C2`, `V2`
- Units must match the supported units (see Supported Units section)
- Both C1 and C2 must use the same unit system (both molar OR both mass)

#### Alternative Format (Legacy: Units in Cells)

You can also use the older format where each cell contains both value and unit:

| C1 | C2 | V2 |
|----|----|-----|
| 10 M | 1 M | 100 mL |
| 50 mg/mL | 5 mg/mL | 200 mL |
| 100 mM | 10 mM | 50 µL |

This format is still supported but the header format is recommended for easier data entry.

#### Running Batch Processing

**Option 1: Using the example script**

```powershell
# Edit batch_process_example.py to set your input filename
python batch_process_example.py
```

**Option 2: Using Python code directly**

```python
from dilution_core import process_excel_dilutions

# Process Excel file with units in headers (recommended)
# Output will be saved as 'dilutions_header_format_results.xlsx'
results = process_excel_dilutions(
    input_file='dilutions_header_format.xlsx',
    output_unit='mL'  # Units for V1 and Dilution columns
)

# Or specify a custom output filename
results = process_excel_dilutions(
    input_file='dilutions_header_format.xlsx',
    output_file='my_results.xlsx',
    output_unit='µL'
)
```

#### Output Format

The program creates a new Excel file with your original columns plus:
- **V1**: Volume of stock solution needed (with unit)
- **Dilution**: Volume of solvent to add (with unit)
- **Error**: Any error messages (empty if successful)

Example output (with header format):

| C1 (M) | C2 (mM) | V2 (mL) | V1 | Dilution | Error |
|--------|---------|---------|----------|----------|-------|
| 10 | 1000 | 100 | 10.00 mL | 90.00 mL | |
| 2 | 500 | 200 | 50.00 mL | 150.00 mL | |
| 0.1 | 10 | 50 | 5.00 mL | 45.00 mL | |

#### Example Excel Templates

**Template 1: Units in Headers (Recommended)**

Create `dilutions_header_format.xlsx`:

| C1 (M) | C2 (mM) | V2 (mL) |
|--------|---------|---------|
| 10 | 1000 | 100 |
| 2 | 500 | 200 |
| 0.1 | 10 | 50 |
| 5 | 1000 | 10 |

**Template 2: Units in Cells (Legacy)**

Create `dilutions_input.xlsx`:

| C1 | C2 | V2 |
|----|----|-----|
| 10 M | 1 M | 100 mL |
| 50 mg/mL | 5 mg/mL | 200 mL |
| 100 mM | 10 mM | 50 µL |
| 2 M | 500 µM | 1 mL |

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
- ✅ Parsing values with units from strings
- ✅ Excel batch processing with valid data
- ✅ Excel error handling and validation

## 📁 Project Structure

```
Dilution-calculator/
├── dilution_core.py           # Core calculation engine with unit conversion & batch processing
├── gui_version.py             # Graphical user interface
├── batch_process_example.py   # Example script for Excel batch processing
├── test_dilution.py           # Comprehensive test suite
├── requirements.txt           # Project dependencies
├── LICENSE                    # License file
└── README.md                  # This file
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

**Parameters:**
- `value` (float): Numerical value to convert
- `from_unit` (str): Source unit (e.g., 'mL', 'µL', 'L')
- `to_unit` (str): Target unit

**Returns:**
- `float`: Converted value

### `convert_concentration(value, from_unit, to_unit)`

Convert concentration from one unit to another.

**Parameters:**
- `value` (float): Numerical value to convert
- `from_unit` (str): Source unit (e.g., 'M', 'mM', 'mg/mL')
- `to_unit` (str): Target unit

**Returns:**
- `float`: Converted value

**Raises:**
- `ValueError`: If conversion between molar and mass units attempted

### `process_excel_dilutions(input_file, output_file=None, output_unit='mL')`

Process batch dilution calculations from an Excel file.

**Parameters:**
- `input_file` (str): Path to input Excel file (.xlsx or .xls)
- `output_file` (str, optional): Path to output Excel file (default: adds '_results' to input filename)
- `output_unit` (str): Unit for output volumes V1 and Dilution (default: 'mL')

**Returns:**
- `pandas.DataFrame`: DataFrame with calculated results including V1, Dilution, and Error columns

**Raises:**
- `ImportError`: If pandas or openpyxl not installed
- `ValueError`: If required columns missing or data format invalid

**Expected Excel Format (Recommended):**
- Column 'C1 (M)' or 'C1_M': Stock concentration values (numeric only)
- Column 'C2 (mM)' or 'C2_mM': Final concentration values (numeric only)
- Column 'V2 (mL)' or 'V2_mL': Final volume values (numeric only)

**Alternative Format (Legacy):**
- Column 'C1': Stock concentration with unit (e.g., "10 M")
- Column 'C2': Final concentration with unit (e.g., "1 mM")
- Column 'V2': Final volume with unit (e.g., "100 mL")

### `validate_dilution(c1, c2, v1, v2, tolerance=0.01)`

Validate that a dilution follows C1V1 = C2V2.

**Parameters:**
- `c1` (float): Stock concentration
- `c2` (float): Final concentration  
- `v1` (float): Stock volume
- `v2` (float): Final volume
- `tolerance` (float): Acceptable relative error (default: 0.01 = 1%)

**Returns:**
- `bool`: True if dilution is valid within tolerance

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

### Issue: "pandas not found" for Excel processing
**Solution**: Install pandas and openpyxl:
```powershell
pip install pandas openpyxl
```

### Issue: Excel file format errors
**Solution**: 
1. Ensure column names are exactly: `C1`, `C2`, `V2` (case-sensitive)
2. Each cell must contain both value and unit with a space: "100 mL" not "100mL"
3. Check that both C1 and C2 use the same unit system (molar OR mass, not mixed)
4. Use supported units only (see Supported Units section)

### Issue: Calculator gives unexpected results
**Solution**: 
1. Check that C1 > C2 (you're diluting, not concentrating)
2. Verify unit selections are correct
3. Ensure input values are numeric

### Issue: Tests failing
**Solution**: Make sure you're in the correct directory and Python can find all modules:
```powershell
python test_dilution.py
```

### Issue: Excel output has ERROR in rows
**Solution**: Check the Error column for specific error messages. Common issues:
- Mixing molar and mass concentration units
- Invalid unit formats
- Missing units in cells
- C1 = 0 (division by zero)

### created with claude sonnet 4.5 in the vscode agent mode
### prompt: 
I want to update the program and add more features. The program is designed to help lab workers perform dilutions easily, so I want to include an option for unit conversion (e.g., mL, µL) for both volume and concentration units.

Another improvement is to add a separate test file that verifies the program’s business logic. Additionally, I plan to create a README file that explains the program, including instructions for installing any dependencies if they exist.
