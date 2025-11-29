"""
PCR Primer Designer GUI - Simple interface for designing PCR primers
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox, filedialog
import threading
from pcr_primer_designer import PCRPrimerOptimizer


class PCRPrimerDesignerGUI:
    """GUI for PCR primer design"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("PCR Primer Designer")
        self.root.geometry("1000x700")
        
        self.optimizer = PCRPrimerOptimizer()
        self.current_results = []
        
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the user interface"""
        
        # Main frame
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(0, weight=1)
        
        # Notebook for tabs
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        # Tab 1: Input
        input_frame = ttk.Frame(notebook, padding="10")
        notebook.add(input_frame, text="Design Primers")
        self.setup_input_tab(input_frame)
        
        # Tab 2: Results
        results_frame = ttk.Frame(notebook, padding="10")
        notebook.add(results_frame, text="Results")
        self.setup_results_tab(results_frame)
        
        # Tab 3: Details
        details_frame = ttk.Frame(notebook, padding="10")
        notebook.add(details_frame, text="Primer Details")
        self.setup_details_tab(details_frame)
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, 
                              relief=tk.SUNKEN, anchor=tk.W)
        status_bar.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        
        # Progress bar - configure style for green color
        style = ttk.Style()
        style.configure("green.Horizontal.TProgressbar", troughcolor='white', 
                       background='green', bordercolor='green', lightcolor='green', darkcolor='green')
        self.progress = ttk.Progressbar(main_frame, mode='indeterminate', style="green.Horizontal.TProgressbar")
        self.progress.grid(row=2, column=0, sticky=(tk.W, tk.E))
    
    def setup_input_tab(self, parent):
        """Setup input tab with scrolling"""
        
        # Configure grid
        parent.grid_rowconfigure(0, weight=1)
        parent.grid_columnconfigure(0, weight=1)
        
        # Create canvas and scrollbar
        canvas = tk.Canvas(parent)
        scrollbar = ttk.Scrollbar(parent, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.grid(row=0, column=0, sticky="nsew")
        scrollbar.grid(row=0, column=1, sticky="ns")
        
        # Mouse wheel scrolling
        def on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", on_mousewheel)
        
        # Template sequence section
        template_frame = ttk.LabelFrame(scrollable_frame, text="Template Sequence", padding="10")
        template_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5, padx=5)
        
        self.template_text = scrolledtext.ScrolledText(template_frame, height=8, width=80, wrap=tk.WORD)
        self.template_text.grid(row=0, column=0, sticky=(tk.W, tk.E))
        
        # Bind text change event to update length
        self.template_text.bind('<KeyRelease>', self.update_template_length)
        
        template_buttons = ttk.Frame(template_frame)
        template_buttons.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        
        # Length display
        self.template_length_var = tk.StringVar(value="Length: 0 bp")
        ttk.Label(template_buttons, textvariable=self.template_length_var, foreground='blue').pack(side=tk.RIGHT, padx=10)
        
        ttk.Button(template_buttons, text="Load from File", 
                  command=self.load_template_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(template_buttons, text="Clear", 
                  command=lambda: self.template_text.delete(1.0, tk.END)).pack(side=tk.LEFT)
        
        # Options section
        options_frame = ttk.LabelFrame(scrollable_frame, text="Primer Design Options", padding="10")
        options_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5, padx=5)
        
        # Forward primer range
        ttk.Label(options_frame, text="Forward Primer Region:", font=('TkDefaultFont', 9, 'bold')).grid(row=0, column=0, padx=5, sticky=tk.W, columnspan=2)
        ttk.Label(options_frame, text="Start (bp):").grid(row=1, column=0, padx=5)
        self.fwd_start_var = tk.StringVar(value="0")
        ttk.Entry(options_frame, textvariable=self.fwd_start_var, width=10).grid(row=1, column=1, padx=5)
        
        ttk.Label(options_frame, text="End (bp):").grid(row=1, column=2, padx=5)
        self.fwd_end_var = tk.StringVar(value="50")
        ttk.Entry(options_frame, textvariable=self.fwd_end_var, width=10).grid(row=1, column=3, padx=5)
        
        ttk.Label(options_frame, text="(leave default to amplify from start)", font=('TkDefaultFont', 8, 'italic')).grid(row=2, column=0, columnspan=4, sticky=tk.W, padx=5)
        
        # Reverse primer range
        ttk.Label(options_frame, text="Reverse Primer Region:", font=('TkDefaultFont', 9, 'bold')).grid(row=3, column=0, padx=5, sticky=tk.W, columnspan=2, pady=(10,0))
        ttk.Label(options_frame, text="Start (bp):").grid(row=4, column=0, padx=5)
        self.rev_start_var = tk.StringVar(value="")
        ttk.Entry(options_frame, textvariable=self.rev_start_var, width=10).grid(row=4, column=1, padx=5)
        
        ttk.Label(options_frame, text="End (bp):").grid(row=4, column=2, padx=5)
        self.rev_end_var = tk.StringVar(value="")
        ttk.Entry(options_frame, textvariable=self.rev_end_var, width=10).grid(row=4, column=3, padx=5)
        
        ttk.Label(options_frame, text="(leave empty to amplify until end)", font=('TkDefaultFont', 8, 'italic')).grid(row=5, column=0, columnspan=4, sticky=tk.W, padx=5)
        
        # Additional parameters
        ttk.Separator(options_frame, orient='horizontal').grid(row=6, column=0, columnspan=4, sticky='ew', pady=10)
        
        # Primer length range
        ttk.Label(options_frame, text="Primer Length:").grid(row=7, column=0, padx=5)
        ttk.Label(options_frame, text="Min:").grid(row=7, column=1, padx=5)
        self.min_len_var = tk.StringVar(value="18")
        ttk.Spinbox(options_frame, from_=15, to=30, textvariable=self.min_len_var, width=8).grid(row=7, column=2, padx=5, sticky=tk.W)
        
        ttk.Label(options_frame, text="Max:").grid(row=7, column=3, padx=5)
        self.max_len_var = tk.StringVar(value="25")
        ttk.Spinbox(options_frame, from_=18, to=35, textvariable=self.max_len_var, width=8).grid(row=7, column=4, padx=5, sticky=tk.W)
        
        # Target Tm
        ttk.Label(options_frame, text="Target Tm:").grid(row=8, column=0, padx=5)
        self.target_tm_var = tk.StringVar(value="60")
        ttk.Spinbox(options_frame, from_=50, to=72, textvariable=self.target_tm_var, width=8).grid(row=8, column=1, padx=5, sticky=tk.W)
        ttk.Label(options_frame, text="°C").grid(row=8, column=2, sticky=tk.W)
        
        # Scan step
        ttk.Label(options_frame, text="Scan Step:").grid(row=9, column=0, padx=5)
        self.step_var = tk.StringVar(value="5")
        ttk.Spinbox(options_frame, from_=1, to=20, textvariable=self.step_var, width=8).grid(row=9, column=1, padx=5, sticky=tk.W)
        ttk.Label(options_frame, text="bp (↑ = faster)").grid(row=9, column=2, sticky=tk.W, padx=5)
        
        # Number of results
        ttk.Label(options_frame, text="Top Results:").grid(row=10, column=0, padx=5)
        self.n_results_var = tk.StringVar(value="10")
        ttk.Spinbox(options_frame, from_=5, to=50, textvariable=self.n_results_var, width=8).grid(row=10, column=1, padx=5, sticky=tk.W)
        
        # Action buttons
        button_frame = ttk.Frame(scrollable_frame)
        button_frame.grid(row=2, column=0, columnspan=2, pady=20, padx=5)
        
        self.design_btn = ttk.Button(button_frame, text="🔍 Design Primers", 
                                     command=self.run_design, style='Accent.TButton')
        self.design_btn.pack(side=tk.LEFT, padx=5)
        
        ttk.Button(button_frame, text="Clear All", 
                  command=self.clear_all).pack(side=tk.LEFT, padx=5)
    
    def setup_results_tab(self, parent):
        """Setup results tab"""
        
        # Tree frame
        tree_frame = ttk.Frame(parent)
        tree_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        
        tree_scroll = ttk.Scrollbar(tree_frame, orient=tk.VERTICAL)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Treeview
        columns = ('Rank', 'Fwd Pos', 'Rev Pos', 'Amplicon', 'Score', 'Fwd Tm', 'Rev Tm', 'Warnings')
        self.results_tree = ttk.Treeview(tree_frame, columns=columns, show='headings', 
                                         yscrollcommand=tree_scroll.set, height=20)
        tree_scroll.config(command=self.results_tree.yview)
        
        # Column headings
        self.results_tree.heading('Rank', text='Rank')
        self.results_tree.heading('Fwd Pos', text='Fwd Pos (bp)')
        self.results_tree.heading('Rev Pos', text='Rev Pos (bp)')
        self.results_tree.heading('Amplicon', text='Amplicon (bp)')
        self.results_tree.heading('Score', text='Score')
        self.results_tree.heading('Fwd Tm', text='Fwd Tm (°C)')
        self.results_tree.heading('Rev Tm', text='Rev Tm (°C)')
        self.results_tree.heading('Warnings', text='Warnings')
        
        # Column widths
        self.results_tree.column('Rank', width=50, anchor=tk.CENTER)
        self.results_tree.column('Fwd Pos', width=100, anchor=tk.CENTER)
        self.results_tree.column('Rev Pos', width=100, anchor=tk.CENTER)
        self.results_tree.column('Amplicon', width=100, anchor=tk.CENTER)
        self.results_tree.column('Score', width=80, anchor=tk.CENTER)
        self.results_tree.column('Fwd Tm', width=100, anchor=tk.CENTER)
        self.results_tree.column('Rev Tm', width=100, anchor=tk.CENTER)
        self.results_tree.column('Warnings', width=80, anchor=tk.CENTER)
        
        self.results_tree.pack(fill=tk.BOTH, expand=True)
        
        # Bind selection
        self.results_tree.bind('<<TreeviewSelect>>', self.on_result_select)
        
        # Buttons
        button_frame = ttk.Frame(parent)
        button_frame.grid(row=1, column=0, pady=10)
        
        ttk.Button(button_frame, text="View Details", 
                  command=self.view_selected_details).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Export Selected", 
                  command=self.export_selected).pack(side=tk.LEFT, padx=5)
    
    def setup_details_tab(self, parent):
        """Setup primer details tab"""
        
        self.details_text = scrolledtext.ScrolledText(parent, height=35, width=100, 
                                                      wrap=tk.WORD, font=('Courier', 9))
        self.details_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        parent.columnconfigure(0, weight=1)
        parent.rowconfigure(0, weight=1)
        
        # Configure tags
        self.details_text.tag_configure('header', font=('Courier', 10, 'bold'))
        self.details_text.tag_configure('subheader', font=('Courier', 9, 'bold'))
        self.details_text.tag_configure('good', foreground='green')
        self.details_text.tag_configure('warning', foreground='orange')
        self.details_text.tag_configure('error', foreground='red')
    
    def load_template_file(self):
        """Load template sequence from file"""
        filename = filedialog.askopenfilename(
            title="Select Template File",
            filetypes=[("All Files", "*.*"), ("FASTA", "*.fasta *.fa"), 
                      ("GenBank", "*.gb *.gbk"), ("Text", "*.txt")]
        )
        if filename:
            try:
                from utils import read_sequence_file
                sequence = read_sequence_file(filename)
                self.template_text.delete(1.0, tk.END)
                self.template_text.insert(1.0, sequence)
                # Update length display
                self.update_template_length()
                self.status_var.set(f"Loaded template: {len(sequence)} bp")
                
                # Auto-set reverse primer defaults
                if not self.rev_start_var.get():
                    self.rev_start_var.set(str(max(0, len(sequence) - 50)))
                if not self.rev_end_var.get():
                    self.rev_end_var.set(str(len(sequence)))
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {e}")
    
    def update_template_length(self, event=None):
        """Update template sequence length display"""
        sequence = self.template_text.get(1.0, tk.END).strip().replace('\n', '').replace(' ', '')
        length = len(sequence)
        self.template_length_var.set(f"Length: {length} bp")
    
    def clear_all(self):
        """Clear all inputs"""
        self.template_text.delete(1.0, tk.END)
        self.results_tree.delete(*self.results_tree.get_children())
        self.details_text.delete(1.0, tk.END)
        self.current_results = []
        self.update_template_length()
        self.status_var.set("Cleared all inputs")
    
    def run_design(self):
        """Run primer design in background thread"""
        
        # Get template
        template_seq = self.template_text.get(1.0, tk.END).strip().replace('\n', '').replace(' ', '').upper()
        
        if not template_seq:
            messagebox.showerror("Error", "Please enter a template sequence")
            return
        
        if len(template_seq) < 50:
            messagebox.showerror("Error", "Template sequence too short (minimum 50 bp)")
            return
        
        # Get parameters
        try:
            fwd_start = int(self.fwd_start_var.get())
            fwd_end = int(self.fwd_end_var.get()) if self.fwd_end_var.get() else None
            rev_start = int(self.rev_start_var.get()) if self.rev_start_var.get() else None
            rev_end = int(self.rev_end_var.get()) if self.rev_end_var.get() else None
            min_len = int(self.min_len_var.get())
            max_len = int(self.max_len_var.get())
            target_tm = float(self.target_tm_var.get())
            step = int(self.step_var.get())
            n_results = int(self.n_results_var.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid parameter values")
            return
        
        # Validate
        if min_len > max_len:
            messagebox.showerror("Error", "Minimum primer length must be ≤ maximum length")
            return
        
        # Disable button
        self.design_btn.config(state='disabled')
        self.progress.start()
        self.status_var.set("Designing primers...")
        
        # Run in thread
        def design_thread():
            try:
                results = self.optimizer.design_pcr_primers(
                    template_seq,
                    fwd_start=fwd_start,
                    fwd_end=fwd_end,
                    rev_start=rev_start,
                    rev_end=rev_end,
                    min_primer_length=min_len,
                    max_primer_length=max_len,
                    target_tm=target_tm,
                    step=step,
                    max_results=n_results
                )
                
                self.root.after(0, lambda: self.display_results(results))
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Error", f"Design failed: {e}"))
                self.root.after(0, lambda: self.status_var.set("Design failed"))
            finally:
                self.root.after(0, lambda: self.progress.stop())
                self.root.after(0, lambda: self.design_btn.config(state='normal'))
        
        thread = threading.Thread(target=design_thread, daemon=True)
        thread.start()
    
    def display_results(self, results):
        """Display design results"""
        
        self.results_tree.delete(*self.results_tree.get_children())
        self.current_results = results
        
        for i, result in enumerate(results, 1):
            self.results_tree.insert('', tk.END, values=(
                i,
                result['fwd_pos'],
                result['rev_pos'],
                result['amplicon_size'],
                f"{result['score']:.1f}",
                f"{result['validations']['forward']['metrics']['tm']:.1f}",
                f"{result['validations']['reverse']['metrics']['tm']:.1f}",
                result['warnings']
            ))
        
        self.status_var.set(f"Found {len(results)} primer pairs")
        
        if results:
            first_item = self.results_tree.get_children()[0]
            self.results_tree.selection_set(first_item)
            self.results_tree.focus(first_item)
    
    def on_result_select(self, event):
        """Handle result selection"""
        selection = self.results_tree.selection()
        if selection:
            item = self.results_tree.item(selection[0])
            rank = int(item['values'][0]) - 1
            if 0 <= rank < len(self.current_results):
                self.display_primer_details(self.current_results[rank])
    
    def display_primer_details(self, result):
        """Display detailed primer information"""
        
        self.details_text.delete(1.0, tk.END)
        
        # Header
        self.details_text.insert(tk.END, "="*80 + "\n", 'header')
        self.details_text.insert(tk.END, "PCR PRIMER DETAILS\n", 'header')
        self.details_text.insert(tk.END, f"Forward Primer Position: {result['fwd_pos']} bp\n", 'header')
        self.details_text.insert(tk.END, f"Reverse Primer Position: {result['rev_pos']} bp\n", 'header')
        self.details_text.insert(tk.END, f"Amplicon Size: {result['amplicon_size']} bp\n", 'header')
        self.details_text.insert(tk.END, f"Overall Score: {result['score']:.1f}/100\n", 'header')
        self.details_text.insert(tk.END, "="*80 + "\n\n", 'header')
        
        primers = result['primers']
        vals = result['validations']
        
        # Forward primer
        self.details_text.insert(tk.END, "FORWARD PRIMER\n", 'subheader')
        self.details_text.insert(tk.END, "-"*80 + "\n")
        self.details_text.insert(tk.END, f"5'- {primers['forward']} -3'\n\n")
        self.display_validation_summary(vals['forward'])
        
        # Reverse primer
        self.details_text.insert(tk.END, "\n\nREVERSE PRIMER\n", 'subheader')
        self.details_text.insert(tk.END, "-"*80 + "\n")
        self.details_text.insert(tk.END, f"5'- {primers['reverse']} -3'\n\n")
        self.display_validation_summary(vals['reverse'])
    
    def display_validation_summary(self, validation):
        """Display validation summary"""
        
        metrics = validation['metrics']
        
        self.details_text.insert(tk.END, f"Length: {metrics['length']} nt\n")
        self.details_text.insert(tk.END, f"GC Content: {metrics['gc_content']}%\n")
        self.details_text.insert(tk.END, f"Tm: {metrics['tm']:.1f}°C\n")
        
        if validation['warnings']:
            self.details_text.insert(tk.END, f"⚠ Warnings: {len(validation['warnings'])}\n", 'warning')
            for warning in validation['warnings']:
                self.details_text.insert(tk.END, f"  • {warning}\n", 'warning')
        else:
            self.details_text.insert(tk.END, "✓ No warnings\n", 'good')
        
        if validation['errors']:
            self.details_text.insert(tk.END, f"✗ Errors: {len(validation['errors'])}\n", 'error')
            for error in validation['errors']:
                self.details_text.insert(tk.END, f"  • {error}\n", 'error')
    
    def view_selected_details(self):
        """Switch to details tab for selected result"""
        selection = self.results_tree.selection()
        if not selection:
            messagebox.showinfo("Info", "Please select a result to view details")
    
    def export_selected(self):
        """Export selected result"""
        selection = self.results_tree.selection()
        if not selection:
            messagebox.showinfo("Info", "Please select a result to export")
            return
        
        item = self.results_tree.item(selection[0])
        rank = int(item['values'][0]) - 1
        
        if 0 <= rank < len(self.current_results):
            self.export_result(self.current_results[rank])
    
    def export_result(self, result):
        """Export primer design result"""
        filename = filedialog.asksaveasfilename(
            title="Export Primer Design",
            defaultextension=".txt",
            filetypes=[("Text", "*.txt"), ("All Files", "*.*")]
        )
        
        if filename:
            try:
                with open(filename, 'w') as f:
                    f.write("PCR Primer Design\n")
                    f.write("="*80 + "\n\n")
                    f.write(f"Forward Primer Position: {result['fwd_pos']} bp\n")
                    f.write(f"Reverse Primer Position: {result['rev_pos']} bp\n")
                    f.write(f"Amplicon Size: {result['amplicon_size']} bp\n")
                    f.write(f"Score: {result['score']:.1f}/100\n\n")
                    
                    f.write("PRIMERS\n")
                    f.write("-"*80 + "\n")
                    f.write(f"Forward: 5'- {result['primers']['forward']} -3'\n")
                    f.write(f"Reverse: 5'- {result['primers']['reverse']} -3'\n\n")
                    
                    f.write("METRICS\n")
                    f.write("-"*80 + "\n")
                    fwd_val = result['validations']['forward']
                    rev_val = result['validations']['reverse']
                    f.write(f"Forward Tm: {fwd_val['metrics']['tm']:.1f}°C\n")
                    f.write(f"Reverse Tm: {rev_val['metrics']['tm']:.1f}°C\n")
                    f.write(f"Forward GC: {fwd_val['metrics']['gc_content']}%\n")
                    f.write(f"Reverse GC: {rev_val['metrics']['gc_content']}%\n")
                
                messagebox.showinfo("Success", f"Exported to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {e}")


def main():
    """Main entry point"""
    root = tk.Tk()
    app = PCRPrimerDesignerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
