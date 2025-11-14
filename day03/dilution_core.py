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