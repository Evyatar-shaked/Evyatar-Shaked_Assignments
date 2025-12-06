"""
Monod Kinetics GUI

GUI application for Monod kinetics calculations and simulations.
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox, filedialog
import sys
from io import StringIO
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from monod_kinetics import MonodKinetics


class MonodKineticsGUI:
    """
    GUI for Monod kinetics calculator.
    """
    
    def __init__(self, root):
        self.root = root
        self.root.title("Monod Kinetics Calculator")
        self.root.geometry("1000x750")
        
        # Configure style
        style = ttk.Style()
        style.configure('Header.TLabel', font=('Arial', 12, 'bold'))
        style.configure('Title.TLabel', font=('Arial', 14, 'bold'))
        
        self.monod = None
        self.results = None
        
        self.create_widgets()
    
    def create_widgets(self):
        """Create all GUI widgets."""
        
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Title
        title = ttk.Label(main_frame, text="Monod Kinetics Calculator", 
                         style='Title.TLabel')
        title.grid(row=0, column=0, columnspan=2, pady=10)
        
        # ===== Parameters Section =====
        params_frame = ttk.LabelFrame(main_frame, text="Monod Parameters", padding="10")
        params_frame.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N), pady=5, padx=5)
        
        # μmax
        ttk.Label(params_frame, text="μmax (hr⁻¹):").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.mu_max_entry = ttk.Entry(params_frame, width=15)
        self.mu_max_entry.insert(0, "0.3")
        self.mu_max_entry.grid(row=0, column=1, sticky=tk.W, pady=5, padx=5)
        ttk.Label(params_frame, text="Maximum specific growth rate", 
                 font=('Arial', 8, 'italic')).grid(row=0, column=2, sticky=tk.W, padx=5)
        
        # Ks
        ttk.Label(params_frame, text="Ks (mM):").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.ks_entry = ttk.Entry(params_frame, width=15)
        self.ks_entry.insert(0, "2.0")
        self.ks_entry.grid(row=1, column=1, sticky=tk.W, pady=5, padx=5)
        ttk.Label(params_frame, text="Half-saturation constant", 
                 font=('Arial', 8, 'italic')).grid(row=1, column=2, sticky=tk.W, padx=5)
        
        # Y_xs
        ttk.Label(params_frame, text="Y_x/s (g/mmol):").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.yxs_entry = ttk.Entry(params_frame, width=15)
        self.yxs_entry.insert(0, "0.02")
        self.yxs_entry.grid(row=2, column=1, sticky=tk.W, pady=5, padx=5)
        ttk.Label(params_frame, text="Yield coefficient (biomass/substrate)", 
                 font=('Arial', 8, 'italic')).grid(row=2, column=2, sticky=tk.W, padx=5)
        
        # Volume
        ttk.Label(params_frame, text="Volume (mL):").grid(row=3, column=0, sticky=tk.W, pady=5)
        self.volume_entry = ttk.Entry(params_frame, width=15)
        self.volume_entry.insert(0, "100")
        self.volume_entry.grid(row=3, column=1, sticky=tk.W, pady=5, padx=5)
        
        # OD conversion
        ttk.Label(params_frame, text="OD conversion (g/L/OD):").grid(row=4, column=0, sticky=tk.W, pady=5)
        self.od_conv_entry = ttk.Entry(params_frame, width=15)
        self.od_conv_entry.insert(0, "0.4")
        self.od_conv_entry.grid(row=4, column=1, sticky=tk.W, pady=5, padx=5)
        ttk.Label(params_frame, text="Biomass per OD unit (typical: 0.4 for E. coli)", 
                 font=('Arial', 8, 'italic')).grid(row=4, column=2, sticky=tk.W, padx=5)
        
        # ===== Simulation Section =====
        sim_frame = ttk.LabelFrame(main_frame, text="Batch Simulation", padding="10")
        sim_frame.grid(row=2, column=0, sticky=(tk.W, tk.E, tk.N), pady=5, padx=5)
        
        # Initial OD
        ttk.Label(sim_frame, text="Initial OD600:").grid(row=0, column=0, sticky=tk.W, pady=5)
        self.init_od_entry = ttk.Entry(sim_frame, width=15)
        self.init_od_entry.insert(0, "0.1")
        self.init_od_entry.grid(row=0, column=1, sticky=tk.W, pady=5, padx=5)
        
        # Initial substrate
        ttk.Label(sim_frame, text="Initial substrate (mM):").grid(row=1, column=0, sticky=tk.W, pady=5)
        self.init_substrate_entry = ttk.Entry(sim_frame, width=15)
        self.init_substrate_entry.insert(0, "20.0")
        self.init_substrate_entry.grid(row=1, column=1, sticky=tk.W, pady=5, padx=5)
        
        # Time
        ttk.Label(sim_frame, text="Simulation time (hours):").grid(row=2, column=0, sticky=tk.W, pady=5)
        self.time_entry = ttk.Entry(sim_frame, width=15)
        self.time_entry.insert(0, "48")
        self.time_entry.grid(row=2, column=1, sticky=tk.W, pady=5, padx=5)
        
        # Buttons
        button_frame = ttk.Frame(sim_frame)
        button_frame.grid(row=3, column=0, columnspan=3, pady=10)
        
        ttk.Button(button_frame, text="Run Simulation", 
                  command=self.run_simulation).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Plot Results", 
                  command=self.plot_results).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Predict Substrate Needed", 
                  command=self.predict_substrate).pack(side=tk.LEFT, padx=5)
        
        # ===== Results Section =====
        results_frame = ttk.LabelFrame(main_frame, text="Results", padding="10")
        results_frame.grid(row=1, column=1, rowspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), 
                          pady=5, padx=5)
        
        self.results_text = scrolledtext.ScrolledText(results_frame, width=60, height=35, 
                                                      font=('Courier', 9))
        self.results_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(1, weight=1)
        main_frame.rowconfigure(2, weight=1)
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)
    
    def get_parameters(self):
        """Get parameters from GUI."""
        try:
            mu_max = float(self.mu_max_entry.get())
            ks = float(self.ks_entry.get())
            yxs = float(self.yxs_entry.get())
            volume = float(self.volume_entry.get())
            
            if mu_max <= 0 or ks <= 0 or yxs <= 0 or volume <= 0:
                messagebox.showerror("Input Error", "All parameters must be positive")
                return None
            
            return mu_max, ks, yxs, volume
        except ValueError:
            messagebox.showerror("Input Error", "Invalid parameter values")
            return None
    
    def run_simulation(self):
        """Run batch simulation."""
        params = self.get_parameters()
        if params is None:
            return
        
        mu_max, ks, yxs, volume = params
        
        try:
            init_od = float(self.init_od_entry.get())
            init_substrate = float(self.init_substrate_entry.get())
            time_hours = float(self.time_entry.get())
            
            if init_od <= 0 or init_substrate < 0 or time_hours <= 0:
                messagebox.showerror("Input Error", "Invalid simulation parameters")
                return
            
            # Create Monod model
            self.monod = MonodKinetics(mu_max, ks, yxs, volume)
            
            # Clear results
            self.results_text.delete(1.0, tk.END)
            
            # Redirect output
            old_stdout = sys.stdout
            sys.stdout = StringIO()
            
            try:
                # Print parameters
                self.monod.print_summary()
                
                # Run simulation
                self.results = self.monod.simulate_batch(init_od, init_substrate, time_hours)
                
                # Print results
                print("\nBATCH SIMULATION RESULTS")
                print("="*70)
                print(f"\nInitial Conditions:")
                print(f"  OD600:              {self.results['OD600'][0]:.3f}")
                print(f"  Substrate:          {self.results['substrate_mM'][0]:.2f} mM")
                print(f"  Biomass:            {self.results['biomass_g_L'][0]:.3f} g/L")
                
                print(f"\nFinal Conditions (after {self.results['time_hr'][-1]:.1f} hours):")
                print(f"  OD600:              {self.results['OD600'][-1]:.3f}")
                print(f"  Substrate:          {self.results['substrate_mM'][-1]:.2f} mM")
                print(f"  Biomass:            {self.results['biomass_g_L'][-1]:.3f} g/L")
                print(f"  Substrate consumed: {self.results['substrate_consumed_mM'][-1]:.2f} mM")
                print(f"  Growth rate:        {self.results['growth_rate_per_hr'][-1]:.4f} hr⁻¹")
                
                # Find maximum OD
                max_od_idx = np.argmax(self.results['OD600'])
                print(f"\nMaximum Growth:")
                print(f"  Time:               {self.results['time_hr'][max_od_idx]:.1f} hours")
                print(f"  OD600:              {self.results['OD600'][max_od_idx]:.3f}")
                print(f"  Fold increase:      {self.results['OD600'][max_od_idx] / self.results['OD600'][0]:.1f}x")
                
                # Doubling time
                print(f"\nGrowth Kinetics:")
                print(f"  Theoretical doubling time at μmax: {np.log(2)/mu_max:.2f} hours")
                if self.results['growth_rate_per_hr'][0] > 0:
                    actual_doubling = np.log(2) / self.results['growth_rate_per_hr'][0]
                    print(f"  Initial doubling time:              {actual_doubling:.2f} hours")
                
                print("\n" + "="*70)
                
                output = sys.stdout.getvalue()
                self.results_text.insert(1.0, output)
                
            finally:
                sys.stdout = old_stdout
            
            messagebox.showinfo("Success", "Simulation completed successfully!")
            
        except Exception as e:
            messagebox.showerror("Error", f"Simulation failed: {str(e)}")
    
    def plot_results(self):
        """Plot simulation results."""
        if self.results is None:
            messagebox.showwarning("No Results", "Please run a simulation first")
            return
        
        # Create plots
        fig, axes = plt.subplots(2, 2, figsize=(12, 8))
        
        time = self.results['time_hr']
        
        # OD vs time
        axes[0, 0].plot(time, self.results['OD600'], 'b-', linewidth=2)
        axes[0, 0].set_xlabel('Time (hours)')
        axes[0, 0].set_ylabel('OD600')
        axes[0, 0].set_title('Bacterial Growth')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Substrate vs time
        axes[0, 1].plot(time, self.results['substrate_mM'], 'g-', linewidth=2)
        axes[0, 1].set_xlabel('Time (hours)')
        axes[0, 1].set_ylabel('Substrate (mM)')
        axes[0, 1].set_title('Substrate Concentration')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Growth rate vs time
        axes[1, 0].plot(time, self.results['growth_rate_per_hr'], 'r-', linewidth=2)
        axes[1, 0].set_xlabel('Time (hours)')
        axes[1, 0].set_ylabel('μ (hr⁻¹)')
        axes[1, 0].set_title('Specific Growth Rate')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Substrate consumed vs time
        axes[1, 1].plot(time, self.results['substrate_consumed_mM'], 'm-', linewidth=2)
        axes[1, 1].set_xlabel('Time (hours)')
        axes[1, 1].set_ylabel('Substrate Consumed (mM)')
        axes[1, 1].set_title('Cumulative Substrate Consumption')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.show()
    
    def predict_substrate(self):
        """Predict substrate needed for target OD."""
        params = self.get_parameters()
        if params is None:
            return
        
        mu_max, ks, yxs, volume = params
        
        # Create dialog for target OD
        dialog = tk.Toplevel(self.root)
        dialog.title("Predict Substrate Needed")
        dialog.geometry("400x250")
        
        ttk.Label(dialog, text="Initial OD600:", font=('Arial', 10)).grid(row=0, column=0, 
                                                                          sticky=tk.W, padx=10, pady=10)
        init_od_entry = ttk.Entry(dialog, width=15)
        init_od_entry.insert(0, "0.1")
        init_od_entry.grid(row=0, column=1, padx=10, pady=10)
        
        ttk.Label(dialog, text="Target OD600:", font=('Arial', 10)).grid(row=1, column=0, 
                                                                         sticky=tk.W, padx=10, pady=10)
        target_od_entry = ttk.Entry(dialog, width=15)
        target_od_entry.insert(0, "1.0")
        target_od_entry.grid(row=1, column=1, padx=10, pady=10)
        
        ttk.Label(dialog, text="Initial substrate (mM):", font=('Arial', 10)).grid(row=2, column=0, 
                                                                                   sticky=tk.W, padx=10, pady=10)
        init_sub_entry = ttk.Entry(dialog, width=15)
        init_sub_entry.insert(0, "10.0")
        init_sub_entry.grid(row=2, column=1, padx=10, pady=10)
        
        def calculate():
            try:
                init_od = float(init_od_entry.get())
                target_od = float(target_od_entry.get())
                init_sub = float(init_sub_entry.get())
                
                monod = MonodKinetics(mu_max, ks, yxs, volume)
                prediction = monod.predict_substrate_needed(init_od, target_od, init_sub)
                
                # Show results
                self.results_text.delete(1.0, tk.END)
                
                output = "\nSUBSTRATE REQUIREMENT PREDICTION\n"
                output += "="*70 + "\n"
                output += f"\nTo grow from OD {prediction['initial_OD']} to OD {prediction['final_OD']}:\n"
                output += f"  Biomass increase:        {prediction['biomass_increase_g_L']:.3f} g/L\n"
                output += f"  Substrate needed:        {prediction['substrate_needed_mM']:.2f} mM\n"
                output += f"  Initial substrate:       {prediction['initial_substrate_mM']:.2f} mM\n"
                output += f"  Final substrate:         {prediction['final_substrate_mM']:.2f} mM\n"
                output += f"  Substrate sufficient:    {prediction['substrate_sufficient']}\n"
                
                if not prediction['substrate_sufficient']:
                    deficit = abs(prediction['final_substrate_mM'])
                    output += f"\n  ⚠ WARNING: Need {deficit:.2f} mM more substrate!\n"
                
                output += "="*70 + "\n"
                
                self.results_text.insert(1.0, output)
                dialog.destroy()
                
            except ValueError:
                messagebox.showerror("Input Error", "Invalid values entered")
        
        ttk.Button(dialog, text="Calculate", command=calculate).grid(row=3, column=0, 
                                                                     columnspan=2, pady=20)


def main():
    """Main entry point."""
    root = tk.Tk()
    app = MonodKineticsGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
