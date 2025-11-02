import tkinter as tk
from tkinter import ttk, messagebox
import math

def circle_area(radius):
    return math.pi * radius ** 2

class CircleAreaGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Circle Area Calculator")
        self.root.geometry("400x250")
        self.root.resizable(False, False)
        
        main_frame = ttk.Frame(root, padding="20")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        ttk.Label(main_frame, text="Circle Area Calculator", 
                  font=("Arial", 16, "bold")).grid(row=0, column=0, columnspan=2, pady=(0, 20))
        
        # Radius input
        ttk.Label(main_frame, text="Radius:", font=("Arial", 12)).grid(row=1, column=0, sticky=tk.W, pady=5)
        self.radius_entry = ttk.Entry(main_frame, font=("Arial", 12), width=15)
        self.radius_entry.grid(row=1, column=1, sticky=(tk.W, tk.E), pady=5, padx=(10, 0))
        
        # Calculate button
        self.calculate_btn = ttk.Button(main_frame, text="Calculate Area", 
                                       command=self.calculate_area)
        self.calculate_btn.grid(row=2, column=0, columnspan=2, pady=20)
        
        # Result display
        ttk.Label(main_frame, text="Result:", font=("Arial", 12, "bold")).grid(row=3, column=0, sticky=tk.W, pady=5)
        self.result_label = ttk.Label(main_frame, text="", font=("Arial", 12), 
                                     foreground="blue", background="lightgray", 
                                     relief="sunken", padding=5)
        self.result_label.grid(row=3, column=1, sticky=(tk.W, tk.E), pady=5, padx=(10, 0))
        
        # Clear button
        self.clear_btn = ttk.Button(main_frame, text="Clear", command=self.clear_fields)
        self.clear_btn.grid(row=4, column=0, columnspan=2, pady=10)
        
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        
        root.bind('<Return>', lambda event: self.calculate_area())
        self.radius_entry.focus()
    
    def calculate_area(self):
        try:
            radius_str = self.radius_entry.get().strip()
            
            if not radius_str:
                messagebox.showerror("Error", "Please enter a radius value.")
                return
            
            radius = float(radius_str)
            if radius <= 0:
                messagebox.showerror("Error", "Radius must be a positive number.")
                return
            
            area = circle_area(radius)
            self.result_label.config(text=f"{area:.2f}")
        
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid number for radius.")
    
    def clear_fields(self):
        self.radius_entry.delete(0, tk.END)
        self.result_label.config(text="")
        self.radius_entry.focus()

def main():
    root = tk.Tk()
    app = CircleAreaGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
