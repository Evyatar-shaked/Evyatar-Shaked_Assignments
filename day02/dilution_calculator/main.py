# main.py
from dilution_core import calculate_dilution

def interactive_mode():
    print("=== Dilution Calculator ===")
    c1 = float(input("Enter stock concentration (C1): "))
    c2 = float(input("Enter final concentration (C2): "))
    v2 = float(input("Enter final volume (V2): "))

    v1, added = calculate_dilution(c1, c2, v2)
    print(f"\nYou need {v1:.3f} units of stock solution.")
    print(f"Add {added:.3f} units of solvent to reach {v2} total volume.")

if __name__ == "__main__":
    interactive_mode()
