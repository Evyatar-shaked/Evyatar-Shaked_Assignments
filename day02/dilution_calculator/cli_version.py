# cli_version.py
import sys
from dilution_core import calculate_dilution

def cli_mode():
    if len(sys.argv) < 4:
        print("Usage: python cli_version.py <C1> <C2> <V2>")
        sys.exit(1)

    c1 = float(sys.argv[1])
    c2 = float(sys.argv[2])
    v2 = float(sys.argv[3])

    v1, added = calculate_dilution(c1, c2, v2)
    print(f"You need {v1:.3f} units of stock solution.")
    print(f"Add {added:.3f} units of solvent to reach {v2} total volume.")

if __name__ == "__main__":
    cli_mode()
