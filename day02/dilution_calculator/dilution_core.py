# dilution_core.py
def calculate_dilution(c1, c2, v2):
    """Return v1 and added volume based on C1V1 = C2V2."""
    if c1 == 0:
        raise ValueError("C1 cannot be zero (division by zero).")

    v1 = (c2 * v2) / c1
    added = v2 - v1
    return v1, added