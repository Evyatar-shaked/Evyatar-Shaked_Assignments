# gui_version.py
"""
GUI version of the dilution calculator with unit conversion support.
Provides an easy-to-use interface for lab workers to calculate dilutions.
"""
import tkinter as tk
from tkinter import ttk, messagebox
from dilution_core import (
    calculate_dilution, 
    VOLUME_UNITS, 
    MOLAR_CONCENTRATION_UNITS,
    MASS_CONCENTRATION_UNITS
)

# Separate unit lists for easier organization
MOLAR_UNITS = list(MOLAR_CONCENTRATION_UNITS.keys())
MASS_UNITS = list(MASS_CONCENTRATION_UNITS.keys())
ALL_CONCENTRATION_UNITS = MOLAR_UNITS + MASS_UNITS


def calculate():
    """Handle the calculate button click event."""
    try:
        c1 = float(entry_c1.get())
        c2 = float(entry_c2.get())
        v2 = float(entry_v2.get())
        
        c1_unit = combo_c1_unit.get()
        c2_unit = combo_c2_unit.get()
        v2_unit = combo_v2_unit.get()
        output_unit = combo_output_unit.get()

        v1, added = calculate_dilution(c1, c2, v2, c1_unit, c2_unit, v2_unit, output_unit)
        
        result_text = (
            f"Stock Volume (V1) = {v1:.4f} {output_unit}\n"
            f"Solvent to Add = {added:.4f} {output_unit}\n\n"
            f"Instructions:\n"
            f"1. Take {v1:.4f} {output_unit} of stock solution (C1)\n"
            f"2. Add {added:.4f} {output_unit} of solvent\n"
            f"3. Final volume will be {v2:.4f} {v2_unit}"
        )
        messagebox.showinfo("Dilution Result", result_text)

    except ValueError as e:
        messagebox.showerror("Error", f"Invalid input: {str(e)}")
    except Exception as e:
        messagebox.showerror("Error", f"Calculation error: {str(e)}")


def clear_fields():
    """Clear all input fields."""
    entry_c1.delete(0, tk.END)
    entry_c2.delete(0, tk.END)
    entry_v2.delete(0, tk.END)


# Create main window
root = tk.Tk()
root.title("Lab Dilution Calculator")
root.geometry("450x500")
root.resizable(False, False)

# Configure styling
style = ttk.Style()
style.configure('Title.TLabel', font=('Arial', 14, 'bold'))
style.configure('Info.TLabel', font=('Arial', 9, 'italic'))

# Title
title_frame = ttk.Frame(root, padding="10")
title_frame.pack(fill=tk.X)
ttk.Label(title_frame, text="🧪 Dilution Calculator", style='Title.TLabel').pack()
ttk.Label(title_frame, text="Based on C1V1 = C2V2", style='Info.TLabel').pack()

# Main frame
main_frame = ttk.Frame(root, padding="20")
main_frame.pack(fill=tk.BOTH, expand=True)

# C1 input
ttk.Label(main_frame, text="Stock Concentration (C1):", font=('Arial', 10, 'bold')).grid(row=0, column=0, sticky=tk.W, pady=5)
c1_frame = ttk.Frame(main_frame)
c1_frame.grid(row=1, column=0, sticky=tk.W, pady=5)
entry_c1 = ttk.Entry(c1_frame, width=15)
entry_c1.pack(side=tk.LEFT, padx=(0, 5))
combo_c1_unit = ttk.Combobox(c1_frame, values=ALL_CONCENTRATION_UNITS, width=10, state='readonly')
combo_c1_unit.set('M')
combo_c1_unit.pack(side=tk.LEFT)

# C2 input
ttk.Label(main_frame, text="Final Concentration (C2):", font=('Arial', 10, 'bold')).grid(row=2, column=0, sticky=tk.W, pady=5)
c2_frame = ttk.Frame(main_frame)
c2_frame.grid(row=3, column=0, sticky=tk.W, pady=5)
entry_c2 = ttk.Entry(c2_frame, width=15)
entry_c2.pack(side=tk.LEFT, padx=(0, 5))
combo_c2_unit = ttk.Combobox(c2_frame, values=ALL_CONCENTRATION_UNITS, width=10, state='readonly')
combo_c2_unit.set('mM')
combo_c2_unit.pack(side=tk.LEFT)

# V2 input
ttk.Label(main_frame, text="Final Volume (V2):", font=('Arial', 10, 'bold')).grid(row=4, column=0, sticky=tk.W, pady=5)
v2_frame = ttk.Frame(main_frame)
v2_frame.grid(row=5, column=0, sticky=tk.W, pady=5)
entry_v2 = ttk.Entry(v2_frame, width=15)
entry_v2.pack(side=tk.LEFT, padx=(0, 5))
combo_v2_unit = ttk.Combobox(v2_frame, values=list(VOLUME_UNITS.keys()), width=8, state='readonly')
combo_v2_unit.set('mL')
combo_v2_unit.pack(side=tk.LEFT)

# Output unit selection
ttk.Label(main_frame, text="Output Volume Unit:", font=('Arial', 10, 'bold')).grid(row=6, column=0, sticky=tk.W, pady=5)
combo_output_unit = ttk.Combobox(main_frame, values=list(VOLUME_UNITS.keys()), width=8, state='readonly')
combo_output_unit.set('mL')
combo_output_unit.grid(row=7, column=0, sticky=tk.W, pady=5)

# Separator
ttk.Separator(main_frame, orient='horizontal').grid(row=8, column=0, sticky='ew', pady=15)

# Buttons
button_frame = ttk.Frame(main_frame)
button_frame.grid(row=9, column=0, pady=10)
ttk.Button(button_frame, text="Calculate", command=calculate, width=15).pack(side=tk.LEFT, padx=5)
ttk.Button(button_frame, text="Clear", command=clear_fields, width=15).pack(side=tk.LEFT, padx=5)

# Info text
info_frame = ttk.Frame(main_frame)
info_frame.grid(row=10, column=0, pady=10)
info_text = (
    "Example:\n"
    "To dilute 10 M stock to 1 M in 100 mL:\n"
    "C1=10, C2=1, V2=100\n"
    "Result: Use 10 mL stock + 90 mL solvent"
)
ttk.Label(info_frame, text=info_text, font=('Arial', 8), justify=tk.LEFT, 
          foreground='gray').pack()

root.mainloop()
