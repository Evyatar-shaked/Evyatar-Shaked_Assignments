# gui_version.py
import tkinter as tk
from tkinter import messagebox
from dilution_core import calculate_dilution

def calculate():
    try:
        c1 = float(entry_c1.get())
        c2 = float(entry_c2.get())
        v2 = float(entry_v2.get())

        v1, added = calculate_dilution(c1, c2, v2)
        messagebox.showinfo("Result", f"V1 = {v1:.3f}\nAdded Volume = {added:.3f}")

    except ValueError:
        messagebox.showerror("Error", "Please enter valid numbers!")

root = tk.Tk()
root.title("Dilution Calculator")

tk.Label(root, text="C1 (stock concentration):").pack()
entry_c1 = tk.Entry(root); entry_c1.pack()

tk.Label(root, text="C2 (final concentration):").pack()
entry_c2 = tk.Entry(root); entry_c2.pack()

tk.Label(root, text="V2 (final volume):").pack()
entry_v2 = tk.Entry(root); entry_v2.pack()

tk.Button(root, text="Calculate", command=calculate).pack(pady=10)

root.mainloop()
