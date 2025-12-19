"""
Formate Beads Calculator - GUI Application

This module provides a graphical user interface for the formate beads calculator.
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import sys
from io import StringIO

# Import core functionality
from formate_beads_core import (
    M07_BEAD_RELEASE, M03_BEAD_RELEASE, 
    ExperimentManager
)


class FormateBeadsGUI:
    """
    GUI Application for Formate Beads Calculator.
    """
    
    def __init__(self, root):
        self.root = root
        self.root.title("Formate Beads Calculator")
        self.root.geometry("900x700")
        
        # Configure style
        style = ttk.Style()
        style.configure('Header.TLabel', font=('Arial', 12, 'bold'))
        style.configure('Title.TLabel', font=('Arial', 14, 'bold'))
        
        self.create_widgets()
    
    def create_widgets(self):
        """Create all GUI widgets."""
        
        # Main container with padding
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Title
        title = ttk.Label(main_frame, text="Formate Beads Release Calculator", 
                         style='Title.TLabel')
        title.grid(row=0, column=0, columnspan=2, pady=10)
        
        # ===== Input Section =====
        input_frame = ttk.LabelFrame(main_frame, text="Experiment Parameters", padding="10")
        input_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        
        # Volume input
        ttk.Label(input_frame, text="Volume (mL):").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.volume_entry = ttk.Entry(input_frame, width=15)
        self.volume_entry.insert(0, "100")
        self.volume_entry.grid(row=0, column=1, sticky=tk.W, pady=5)
        
        # Info label
        info_label = ttk.Label(input_frame, text="Strategy: Mixed M07 + M03 Beads (Cumulative Mode)", 
                              font=('Arial', 9, 'italic'))
        info_label.grid(row=1, column=0, columnspan=2, sticky=tk.W, pady=5)
        
        # ===== Concentration Input Section =====
        conc_frame = ttk.LabelFrame(main_frame, text="Desired Cumulative Concentrations (mM)", padding="10")
        conc_frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        
        ttk.Label(conc_frame, text="Enter desired CUMULATIVE formate concentration for each day:", 
                 font=('Arial', 9, 'italic')).grid(row=0, column=0, columnspan=4, pady=5)
        
        # Create entries for each day
        self.day_entries = {}
        for i in range(7):
            day = i + 1
            col = i % 4
            row = 1 + (i // 4)
            
            ttk.Label(conc_frame, text=f"Day {day}:").grid(row=row, column=col*2, 
                                                           sticky=tk.W, padx=5, pady=3)
            entry = ttk.Entry(conc_frame, width=10)
            entry.insert(0, str(day * 5.0))  # Example cumulative values
            entry.grid(row=row, column=col*2+1, sticky=tk.W, padx=5, pady=3)
            self.day_entries[day] = entry
        
        # ===== Release Profile Configuration =====
        profile_frame = ttk.LabelFrame(main_frame, text="Bead Release Profiles (mmol/day)", 
                                      padding="10")
        profile_frame.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        
        # M07 beads profile
        ttk.Label(profile_frame, text="M07 Beads:", font=('Arial', 9, 'bold')).grid(
            row=0, column=0, sticky=tk.W, pady=2)
        m07_frame = ttk.Frame(profile_frame)
        m07_frame.grid(row=1, column=0, sticky=tk.W)
        
        self.m07_entries = {}
        for day in range(1, 8):
            ttk.Label(m07_frame, text=f"D{day}:").grid(row=0, column=(day-1)*2, 
                                                         sticky=tk.W, padx=2)
            entry = ttk.Entry(m07_frame, width=7)
            entry.insert(0, str(M07_BEAD_RELEASE[day]))
            entry.grid(row=0, column=(day-1)*2+1, padx=2)
            self.m07_entries[day] = entry
        
        # M03 beads profile
        ttk.Label(profile_frame, text="M03 Beads:", font=('Arial', 9, 'bold')).grid(
            row=2, column=0, sticky=tk.W, pady=(10, 2))
        m03_frame = ttk.Frame(profile_frame)
        m03_frame.grid(row=3, column=0, sticky=tk.W)
        
        self.m03_entries = {}
        for day in range(1, 8):
            ttk.Label(m03_frame, text=f"D{day}:").grid(row=0, column=(day-1)*2, 
                                                         sticky=tk.W, padx=2)
            entry = ttk.Entry(m03_frame, width=7)
            entry.insert(0, str(M03_BEAD_RELEASE[day]))
            entry.grid(row=0, column=(day-1)*2+1, padx=2)
            self.m03_entries[day] = entry
        
        # ===== Calculate Button =====
        calc_button = ttk.Button(main_frame, text="Calculate Bead Schedule", 
                                command=self.calculate)
        calc_button.grid(row=4, column=0, columnspan=2, pady=15)
        
        # ===== Results Section =====
        results_frame = ttk.LabelFrame(main_frame, text="Results", padding="10")
        results_frame.grid(row=5, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), 
                          pady=5)
        
        self.results_text = scrolledtext.ScrolledText(results_frame, width=100, height=20, 
                                                      font=('Courier', 9))
        self.results_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights for resizing
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(5, weight=1)
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)
    
    def get_release_profile(self, bead_type):
        """Get release profile from GUI entries."""
        if bead_type == 'M07':
            entries = self.m07_entries
        else:
            entries = self.m03_entries
        
        profile = {}
        for day, entry in entries.items():
            try:
                profile[day] = float(entry.get())
            except ValueError:
                messagebox.showerror("Input Error", 
                                   f"Invalid release rate for {bead_type} bead day {day}")
                return None
        return profile
    
    def calculate(self):
        """Calculate bead schedule based on user inputs."""
        try:
            # Get volume
            volume = float(self.volume_entry.get())
            if volume <= 0:
                messagebox.showerror("Input Error", "Volume must be positive")
                return
            
            # Get desired concentrations
            desired_conc = {}
            for day, entry in self.day_entries.items():
                try:
                    conc = float(entry.get())
                    if conc < 0:
                        messagebox.showerror("Input Error", 
                                           f"Concentration for day {day} cannot be negative")
                        return
                    desired_conc[day] = conc
                except ValueError:
                    messagebox.showerror("Input Error", 
                                       f"Invalid concentration for day {day}")
                    return
            
            # Get release profiles
            m07_profile = self.get_release_profile('M07')
            m03_profile = self.get_release_profile('M03')
            
            if m07_profile is None or m03_profile is None:
                return
            
            # Update global profiles
            import formate_beads_core
            formate_beads_core.M07_BEAD_RELEASE = m07_profile
            formate_beads_core.M03_BEAD_RELEASE = m03_profile
            
            # Clear results
            self.results_text.delete(1.0, tk.END)
            
            # Redirect print to results text widget
            old_stdout = sys.stdout
            sys.stdout = StringIO()
            
            try:
                # Create experiment and calculate
                experiment = ExperimentManager(volume_ml=volume)
                schedule = experiment.calculate_beads_needed(desired_conc)
                
                # Get printed output
                output = sys.stdout.getvalue()
                
                # Display results
                self.results_text.insert(1.0, output)
                
            finally:
                sys.stdout = old_stdout
            
        except ValueError as e:
            messagebox.showerror("Input Error", f"Invalid input: {str(e)}")
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")


def main():
    """Main entry point for GUI application."""
    root = tk.Tk()
    app = FormateBeadsGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
