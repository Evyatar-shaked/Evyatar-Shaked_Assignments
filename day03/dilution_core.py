# dilution_core.py
"""
Core dilution calculation module with unit conversion support.
Provides functions for C1V1 = C2V2 calculations with automatic unit handling.
"""

# Unit conversion factors (to base unit: liters for volume, M for concentration)
VOLUME_UNITS = {
    'L': 1.0,
    'mL': 1e-3,
    'µL': 1e-6,
    'uL': 1e-6,  # Alternative spelling for microliter
    'nL': 1e-9
}

# Molar concentration units (interconvertible among themselves)
MOLAR_CONCENTRATION_UNITS = {
    'M': 1.0,
    'mM': 1e-3,
    'µM': 1e-6,
    'uM': 1e-6,  # Alternative spelling for micromolar
    'nM': 1e-9
}

# Mass concentration units (interconvertible among themselves)
MASS_CONCENTRATION_UNITS = {
    'g/mL': 1.0,
    'mg/mL': 1e-3,
    'µg/mL': 1e-6,
    'ug/mL': 1e-6,  # Alternative spelling for microgram
    'ng/mL': 1e-9,
    'g/L': 1e-3,  # 1 g/L = 0.001 g/mL
    'mg/L': 1e-6,
    'µg/L': 1e-9,
    'ug/L': 1e-9
}

# Combined dictionary for backward compatibility
CONCENTRATION_UNITS = {**MOLAR_CONCENTRATION_UNITS, **MASS_CONCENTRATION_UNITS}


def convert_volume(value, from_unit, to_unit):
    """
    Convert volume from one unit to another.
    
    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'mL', 'µL', 'L')
        to_unit: Target unit
        
    Returns:
        Converted value
        
    Raises:
        ValueError: If units are not recognized
    """
    if from_unit not in VOLUME_UNITS:
        raise ValueError(f"Unknown volume unit: {from_unit}")
    if to_unit not in VOLUME_UNITS:
        raise ValueError(f"Unknown volume unit: {to_unit}")
    
    # Convert to base unit (L) then to target unit
    base_value = value * VOLUME_UNITS[from_unit]
    return base_value / VOLUME_UNITS[to_unit]


def convert_concentration(value, from_unit, to_unit, molecular_weight=None):
    """
    Convert concentration from one unit to another.
    
    Args:
        value: Numerical value to convert
        from_unit: Source unit (e.g., 'M', 'mM', 'µM', 'mg/mL', 'µg/mL')
        to_unit: Target unit
        molecular_weight: Required for conversions between molar and mass units (g/mol)
        
    Returns:
        Converted value
        
    Raises:
        ValueError: If units are not recognized or conversion not possible
    """
    if from_unit not in CONCENTRATION_UNITS:
        raise ValueError(f"Unknown concentration unit: {from_unit}")
    if to_unit not in CONCENTRATION_UNITS:
        raise ValueError(f"Unknown concentration unit: {to_unit}")
    
    # Check if both units are in the same system
    from_is_molar = from_unit in MOLAR_CONCENTRATION_UNITS
    to_is_molar = to_unit in MOLAR_CONCENTRATION_UNITS
    from_is_mass = from_unit in MASS_CONCENTRATION_UNITS
    to_is_mass = to_unit in MASS_CONCENTRATION_UNITS
    
    # Molar to molar conversion
    if from_is_molar and to_is_molar:
        base_value = value * MOLAR_CONCENTRATION_UNITS[from_unit]
        return base_value / MOLAR_CONCENTRATION_UNITS[to_unit]
    
    # Mass to mass conversion
    if from_is_mass and to_is_mass:
        base_value = value * MASS_CONCENTRATION_UNITS[from_unit]
        return base_value / MASS_CONCENTRATION_UNITS[to_unit]
    
    # Cross-system conversions require molecular weight
    raise ValueError(f"Cannot convert between molar ({from_unit}) and mass ({to_unit}) units without molecular weight")


def calculate_dilution(c1, c2, v2, c1_unit='M', c2_unit='M', v2_unit='mL', output_unit='mL'):
    """
    Calculate dilution volumes based on C1V1 = C2V2 with unit support.
    
    Args:
        c1: Stock concentration
        c2: Final concentration
        v2: Final volume
        c1_unit: Unit for c1 (default: 'M')
        c2_unit: Unit for c2 (default: 'M')
        v2_unit: Unit for v2 (default: 'mL')
        output_unit: Desired unit for output volumes (default: 'mL')
        
    Returns:
        tuple: (v1, added_volume) both in output_unit
        
    Raises:
        ValueError: If c1 is zero, units are invalid, or mixing molar and mass units
    """
    if c1 == 0:
        raise ValueError("C1 cannot be zero (division by zero).")
    
    # Check that both concentrations are in the same system (molar or mass)
    c1_is_molar = c1_unit in MOLAR_CONCENTRATION_UNITS
    c2_is_molar = c2_unit in MOLAR_CONCENTRATION_UNITS
    c1_is_mass = c1_unit in MASS_CONCENTRATION_UNITS
    c2_is_mass = c2_unit in MASS_CONCENTRATION_UNITS
    
    if (c1_is_molar and c2_is_mass) or (c1_is_mass and c2_is_molar):
        raise ValueError(f"Cannot mix molar units ({c1_unit}) and mass units ({c2_unit}). Both concentrations must use the same unit system.")
    
    # Convert concentrations to normalized values
    if c1_is_molar and c2_is_molar:
        c1_normalized = c1 * MOLAR_CONCENTRATION_UNITS[c1_unit]
        c2_normalized = c2 * MOLAR_CONCENTRATION_UNITS[c2_unit]
    elif c1_is_mass and c2_is_mass:
        c1_normalized = c1 * MASS_CONCENTRATION_UNITS[c1_unit]
        c2_normalized = c2 * MASS_CONCENTRATION_UNITS[c2_unit]
    else:
        raise ValueError(f"Unknown concentration units: {c1_unit}, {c2_unit}")
    
    # Convert v2 to base unit (L)
    v2_normalized = v2 * VOLUME_UNITS.get(v2_unit, 1.0)
    
    # Calculate v1 in base unit (L)
    v1_normalized = (c2_normalized * v2_normalized) / c1_normalized
    
    # Convert results to desired output unit
    v1 = v1_normalized / VOLUME_UNITS.get(output_unit, 1.0)
    v2_output = v2_normalized / VOLUME_UNITS.get(output_unit, 1.0)
    added = v2_output - v1
    
    return v1, added


def validate_dilution(c1, c2, v1, v2, tolerance=0.01):
    """
    Validate that a dilution follows C1V1 = C2V2.
    
    Args:
        c1: Stock concentration
        c2: Final concentration
        v1: Stock volume
        v2: Final volume
        tolerance: Acceptable relative error (default: 0.01 = 1%)
        
    Returns:
        bool: True if dilution is valid within tolerance
    """
    left_side = c1 * v1
    right_side = c2 * v2
    
    if right_side == 0:
        return left_side == 0
    
    relative_error = abs(left_side - right_side) / right_side
    return relative_error <= tolerance


def parse_value_with_unit(cell_value):
    """
    Parse a cell value that contains both a number and unit.
    
    Args:
        cell_value: String or numeric value (e.g., "100 mL", "10 M", or just 100)
        
    Returns:
        tuple: (value, unit) where value is float and unit is string
        
    Raises:
        ValueError: If the format is invalid
    """
    if isinstance(cell_value, (int, float)):
        raise ValueError(f"Value must include a unit. Got: {cell_value}")
    
    cell_value = str(cell_value).strip()
    parts = cell_value.split()
    
    if len(parts) != 2:
        raise ValueError(f"Invalid format. Expected 'value unit' (e.g., '100 mL'). Got: {cell_value}")
    
    try:
        value = float(parts[0])
        unit = parts[1]
        return value, unit
    except ValueError:
        raise ValueError(f"Could not parse numeric value from: {cell_value}")


def parse_column_header(header):
    """
    Parse a column header to extract the base name and unit.
    
    Args:
        header: Column header string (e.g., "C1 (M)", "C1_M", "C1-M", or just "C1")
        
    Returns:
        tuple: (base_name, unit) where unit is None if not specified
        
    Examples:
        "C1 (M)" -> ("C1", "M")
        "C1_M" -> ("C1", "M")
        "V2-mL" -> ("V2", "mL")
        "C1" -> ("C1", None)
    """
    import re
    
    header = str(header).strip()
    
    # Try different patterns: "C1 (M)", "C1_M", "C1-M", "C1 M"
    # Pattern with parentheses: C1 (M) or C1(M)
    match = re.match(r'^([A-Za-z0-9]+)\s*\(([^)]+)\)$', header)
    if match:
        return match.group(1), match.group(2)
    
    # Pattern with underscore: C1_M
    match = re.match(r'^([A-Za-z0-9]+)_(.+)$', header)
    if match:
        return match.group(1), match.group(2)
    
    # Pattern with dash: C1-M
    match = re.match(r'^([A-Za-z0-9]+)-(.+)$', header)
    if match:
        return match.group(1), match.group(2)
    
    # Pattern with space: C1 M (but only if second part looks like a unit)
    match = re.match(r'^([A-Za-z0-9]+)\s+([A-Za-z/µμ]+)$', header)
    if match:
        return match.group(1), match.group(2)
    
    # No unit specified
    return header, None


def process_excel_dilutions(input_file, output_file=None, output_unit='mL'):
    """
    Process batch dilution calculations from an Excel file.
    
    Expected Excel format (units in column headers):
    - Column 'C1 (M)' or 'C1_M' or 'C1-M': Stock concentration values
    - Column 'C2 (mM)' or 'C2_mM': Final concentration values
    - Column 'V2 (mL)' or 'V2_mL': Final volume values
    
    Alternative format (units in cells - legacy):
    - Column 'C1': Stock concentration with unit (e.g., "10 M")
    - Column 'C2': Final concentration with unit (e.g., "1 mM")
    - Column 'V2': Final volume with unit (e.g., "100 mL")
    
    The function will add two columns:
    - Column 'V1': Required stock volume
    - Column 'Dilution': Volume of solvent to add
    
    Args:
        input_file: Path to input Excel file (.xlsx or .xls)
        output_file: Path to output Excel file (default: adds '_results' to input filename)
        output_unit: Unit for output volumes V1 and Dilution (default: 'mL')
        
    Returns:
        pandas.DataFrame: DataFrame with calculated results
        
    Raises:
        ImportError: If pandas or openpyxl not installed
        ValueError: If required columns missing or data format invalid
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for Excel processing. Install with: pip install pandas openpyxl")
    
    # Read the Excel file
    try:
        df = pd.read_excel(input_file)
    except Exception as e:
        raise ValueError(f"Could not read Excel file: {e}")
    
    # Parse column headers to find C1, C2, V2 and their units
    c1_col = None
    c2_col = None
    v2_col = None
    c1_unit = None
    c2_unit = None
    v2_unit = None
    
    for col in df.columns:
        base_name, unit = parse_column_header(col)
        if base_name == 'C1':
            c1_col = col
            c1_unit = unit
        elif base_name == 'C2':
            c2_col = col
            c2_unit = unit
        elif base_name == 'V2':
            v2_col = col
            v2_unit = unit
    
    # Validate required columns exist
    missing = []
    if c1_col is None:
        missing.append('C1')
    if c2_col is None:
        missing.append('C2')
    if v2_col is None:
        missing.append('V2')
    
    if missing:
        raise ValueError(f"Missing required columns: {missing}. Excel must have columns with base names: C1, C2, V2 (e.g., 'C1 (M)', 'C1_M', or just 'C1')")
    
    # Initialize result columns
    v1_values = []
    dilution_values = []
    errors = []
    
    # Process each row
    for idx, row in df.iterrows():
        try:
            # Get values from cells
            c1_val = row[c1_col]
            c2_val = row[c2_col]
            v2_val = row[v2_col]
            
            # Determine units (either from header or from cell value)
            if c1_unit is None:
                # Legacy format: unit in cell value
                c1_val, c1_unit_row = parse_value_with_unit(c1_val)
            else:
                # New format: unit in header, only numeric value in cell
                c1_val = float(c1_val)
                c1_unit_row = c1_unit
            
            if c2_unit is None:
                c2_val, c2_unit_row = parse_value_with_unit(c2_val)
            else:
                c2_val = float(c2_val)
                c2_unit_row = c2_unit
            
            if v2_unit is None:
                v2_val, v2_unit_row = parse_value_with_unit(v2_val)
            else:
                v2_val = float(v2_val)
                v2_unit_row = v2_unit
            
            # Calculate dilution
            v1, added = calculate_dilution(
                c1=c1_val,
                c2=c2_val,
                v2=v2_val,
                c1_unit=c1_unit_row,
                c2_unit=c2_unit_row,
                v2_unit=v2_unit_row,
                output_unit=output_unit
            )
            
            v1_values.append(f"{v1:.2f} {output_unit}")
            dilution_values.append(f"{added:.2f} {output_unit}")
            errors.append(None)
            
        except Exception as e:
            v1_values.append("ERROR")
            dilution_values.append("ERROR")
            errors.append(str(e))
    
    # Add results to dataframe
    df['V1'] = v1_values
    df['Dilution'] = dilution_values
    df['Error'] = errors
    
    # Save to output file if specified
    if output_file:
        df.to_excel(output_file, index=False)
        print(f"Results saved to: {output_file}")
    elif output_file is None:
        # Generate default output filename
        import os
        base_name = os.path.splitext(input_file)[0]
        output_file = f"{base_name}_results.xlsx"
        df.to_excel(output_file, index=False)
        print(f"Results saved to: {output_file}")
    
    return df