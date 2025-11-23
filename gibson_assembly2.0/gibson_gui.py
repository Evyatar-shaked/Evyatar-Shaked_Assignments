"""
Gibson Assembly Primer Designer - GUI Application
Graphical interface with automatic primer optimization
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox, filedialog
from optimizer import PrimerOptimizer
from primer_design import PrimerDesigner
import threading


class GibsonAssemblyGUI:
    """GUI for Gibson Assembly Primer Designer with optimization"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("Gibson Assembly Primer Designer")
        self.root.geometry("1200x800")
        
        self.optimizer = PrimerOptimizer()
        self.designer = PrimerDesigner()
        self.current_results = []
        
        # Store original sequences and file paths for GenBank export
        self.vector_seq = None
        self.insert_seq = None
        self.vector_file = None
        self.insert_file = None
        
        self.setup_ui()
    
    def setup_ui(self):
        """Setup the user interface"""
        
        # Main container with padding
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(3, weight=1)
        
        # Title
        title_label = ttk.Label(main_frame, text="🧬 Gibson Assembly Primer Designer", 
                               font=('Arial', 16, 'bold'))
        title_label.grid(row=0, column=0, pady=10)
        
        # Create notebook for tabs
        notebook = ttk.Notebook(main_frame)
        notebook.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        # Tab 1: Input
        input_frame = ttk.Frame(notebook, padding="10")
        notebook.add(input_frame, text="Input Sequences")
        self.setup_input_tab(input_frame)
        
        # Tab 2: Results
        results_frame = ttk.Frame(notebook, padding="10")
        notebook.add(results_frame, text="Results")
        self.setup_results_tab(results_frame)
        
        # Tab 3: Primer Details
        details_frame = ttk.Frame(notebook, padding="10")
        notebook.add(details_frame, text="Primer Details")
        self.setup_details_tab(details_frame)
        
        # Status bar
        self.status_var = tk.StringVar(value="Ready")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, 
                              relief=tk.SUNKEN, anchor=tk.W)
        status_bar.grid(row=2, column=0, sticky=(tk.W, tk.E), pady=5)
    
    def setup_input_tab(self, parent):
        """Setup input tab with scrollbar"""
        
        # Configure parent grid
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
        
        # Enable mouse wheel scrolling
        def on_mousewheel(event):
            canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        canvas.bind_all("<MouseWheel>", on_mousewheel)
        
        # Vector section
        vector_frame = ttk.LabelFrame(scrollable_frame, text="Vector Sequence", padding="10")
        vector_frame.grid(row=0, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5, padx=5)
        
        self.vector_text = scrolledtext.ScrolledText(vector_frame, height=8, width=80, wrap=tk.WORD)
        self.vector_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Bind text change event to update length
        self.vector_text.bind('<KeyRelease>', self.update_vector_length)
        
        vector_buttons = ttk.Frame(vector_frame)
        vector_buttons.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        
        # Length display
        self.vector_length_var = tk.StringVar(value="Length: 0 bp")
        ttk.Label(vector_buttons, textvariable=self.vector_length_var, foreground='blue').pack(side=tk.RIGHT, padx=10)
        
        ttk.Button(vector_buttons, text="Load from File", 
                  command=self.load_vector_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(vector_buttons, text="Clear", 
                  command=lambda: self.vector_text.delete(1.0, tk.END)).pack(side=tk.LEFT)
        
        # Insert section
        insert_frame = ttk.LabelFrame(scrollable_frame, text="Insert Sequence", padding="10")
        insert_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5, padx=5)
        
        self.insert_text = scrolledtext.ScrolledText(insert_frame, height=8, width=80, wrap=tk.WORD)
        self.insert_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Bind text change event to update length
        self.insert_text.bind('<KeyRelease>', self.update_insert_length)
        
        insert_buttons = ttk.Frame(insert_frame)
        insert_buttons.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        
        # Length display
        self.insert_length_var = tk.StringVar(value="Length: 0 bp")
        ttk.Label(insert_buttons, textvariable=self.insert_length_var, foreground='blue').pack(side=tk.RIGHT, padx=10)
        
        ttk.Button(insert_buttons, text="Load from File", 
                  command=self.load_insert_file).pack(side=tk.LEFT, padx=5)
        ttk.Button(insert_buttons, text="Clear", 
                  command=lambda: self.insert_text.delete(1.0, tk.END)).pack(side=tk.LEFT)
        
        # Options section
        options_frame = ttk.LabelFrame(scrollable_frame, text="Options", padding="10")
        options_frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5, padx=5)
        
        # Range selection
        range_frame = ttk.Frame(options_frame)
        range_frame.grid(row=0, column=0, columnspan=4, sticky=(tk.W, tk.E), pady=5)
        
        # Vector upstream primer range
        ttk.Label(range_frame, text="Vector Forward Primer Region:", font=('TkDefaultFont', 9, 'bold')).grid(row=0, column=0, padx=5, sticky=tk.W, columnspan=2)
        ttk.Label(range_frame, text="Start (bp):").grid(row=1, column=0, padx=5)
        self.vector_upstream_start_var = tk.StringVar(value="100")
        ttk.Entry(range_frame, textvariable=self.vector_upstream_start_var, width=10).grid(row=1, column=1, padx=5)
        
        ttk.Label(range_frame, text="End (bp):").grid(row=1, column=2, padx=5)
        self.vector_upstream_end_var = tk.StringVar(value="300")
        ttk.Entry(range_frame, textvariable=self.vector_upstream_end_var, width=10).grid(row=1, column=3, padx=5)
        
        # Vector downstream primer range
        ttk.Label(range_frame, text="Vector Reverse Primer Region:", font=('TkDefaultFont', 9, 'bold')).grid(row=2, column=0, padx=5, sticky=tk.W, columnspan=2, pady=(10,0))
        ttk.Label(range_frame, text="Start (bp):").grid(row=3, column=0, padx=5)
        self.vector_downstream_start_var = tk.StringVar(value="700")
        ttk.Entry(range_frame, textvariable=self.vector_downstream_start_var, width=10).grid(row=3, column=1, padx=5)
        
        ttk.Label(range_frame, text="End (bp):").grid(row=3, column=2, padx=5)
        self.vector_downstream_end_var = tk.StringVar(value="900")
        ttk.Entry(range_frame, textvariable=self.vector_downstream_end_var, width=10).grid(row=3, column=3, padx=5)
        
        # Insert optimization toggle
        ttk.Separator(range_frame, orient='horizontal').grid(row=4, column=0, columnspan=4, sticky='ew', pady=10)
        self.optimize_insert_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(range_frame, text="Optimize Insert Primer Regions", 
                       variable=self.optimize_insert_var,
                       command=self.toggle_insert_range).grid(row=5, column=0, columnspan=2, sticky=tk.W, padx=5)
        
        # Insert forward primer range
        ttk.Label(range_frame, text="Insert Forward Primer Region:").grid(row=6, column=0, padx=5, sticky=tk.W, columnspan=2)
        ttk.Label(range_frame, text="Start (bp):").grid(row=7, column=0, padx=5)
        self.insert_fwd_start_var = tk.StringVar(value="0")
        self.insert_fwd_start_entry = ttk.Entry(range_frame, textvariable=self.insert_fwd_start_var, 
                                           width=10, state='disabled')
        self.insert_fwd_start_entry.grid(row=7, column=1, padx=5)
        
        ttk.Label(range_frame, text="End (bp):").grid(row=7, column=2, padx=5)
        self.insert_fwd_end_var = tk.StringVar(value="100")
        self.insert_fwd_end_entry = ttk.Entry(range_frame, textvariable=self.insert_fwd_end_var, 
                                         width=10, state='disabled')
        self.insert_fwd_end_entry.grid(row=7, column=3, padx=5)
        
        # Insert reverse primer range
        ttk.Label(range_frame, text="Insert Reverse Primer Region:").grid(row=8, column=0, padx=5, sticky=tk.W, columnspan=2)
        ttk.Label(range_frame, text="Start (bp):").grid(row=9, column=0, padx=5)
        self.insert_rev_start_var = tk.StringVar(value="900")
        self.insert_rev_start_entry = ttk.Entry(range_frame, textvariable=self.insert_rev_start_var, 
                                           width=10, state='disabled')
        self.insert_rev_start_entry.grid(row=9, column=1, padx=5)
        
        ttk.Label(range_frame, text="End (bp):").grid(row=9, column=2, padx=5)
        self.insert_rev_end_var = tk.StringVar(value="1000")
        self.insert_rev_end_entry = ttk.Entry(range_frame, textvariable=self.insert_rev_end_var, 
                                         width=10, state='disabled')
        self.insert_rev_end_entry.grid(row=9, column=3, padx=5)
        
        # Additional parameters
        ttk.Separator(range_frame, orient='horizontal').grid(row=10, column=0, columnspan=4, sticky='ew', pady=10)
        
        # Homology length
        ttk.Label(range_frame, text="Homology Length:").grid(row=11, column=0, padx=5)
        self.homology_var = tk.StringVar(value="25")
        ttk.Spinbox(range_frame, from_=15, to=40, textvariable=self.homology_var, 
                   width=10).grid(row=11, column=1, padx=5, sticky=tk.W)
        ttk.Label(range_frame, text="nt").grid(row=11, column=2, sticky=tk.W)
        
        # Scan step
        ttk.Label(range_frame, text="Scan Step:").grid(row=12, column=0, padx=5)
        self.step_var = tk.StringVar(value="20")
        ttk.Spinbox(range_frame, from_=5, to=50, textvariable=self.step_var, 
                   width=10).grid(row=12, column=1, padx=5, sticky=tk.W)
        ttk.Label(range_frame, text="bp (↑ = faster)").grid(row=12, column=2, sticky=tk.W, padx=5)
        
        # Number of results
        ttk.Label(range_frame, text="Top Results:").grid(row=13, column=0, padx=5)
        self.n_results_var = tk.StringVar(value="5")
        ttk.Spinbox(range_frame, from_=1, to=20, textvariable=self.n_results_var, 
                   width=10).grid(row=13, column=1, padx=5, sticky=tk.W)
        
        # Action buttons
        button_frame = ttk.Frame(scrollable_frame)
        button_frame.grid(row=3, column=0, columnspan=2, pady=20, padx=5)
        
        self.optimize_btn = ttk.Button(button_frame, text="🔍 Find Best Cut Sites", 
                                       command=self.run_optimization, style='Accent.TButton')
        self.optimize_btn.pack(side=tk.LEFT, padx=10)
        
        ttk.Button(button_frame, text="Clear All", 
                  command=self.clear_all).pack(side=tk.LEFT, padx=10)
        
        # Progress bar
        self.progress = ttk.Progressbar(parent, mode='indeterminate')
        self.progress.grid(row=4, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
    
    def setup_results_tab(self, parent):
        """Setup results tab"""
        
        # Results tree
        tree_frame = ttk.Frame(parent)
        tree_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        # Scrollbar
        tree_scroll = ttk.Scrollbar(tree_frame)
        tree_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        # Treeview
        columns = ('Rank', 'Vector Fwd', 'Vector Rev', 'Insert Fwd', 'Insert Rev', 'Score', 'Warnings', 'Errors')
        self.results_tree = ttk.Treeview(tree_frame, columns=columns, show='headings', 
                                         yscrollcommand=tree_scroll.set, height=20)
        tree_scroll.config(command=self.results_tree.yview)
        
        # Column headings
        self.results_tree.heading('Rank', text='Rank')
        self.results_tree.heading('Vector Fwd', text='Vector Fwd (bp)')
        self.results_tree.heading('Vector Rev', text='Vector Rev (bp)')
        self.results_tree.heading('Insert Fwd', text='Insert Fwd (bp)')
        self.results_tree.heading('Insert Rev', text='Insert Rev (bp)')
        self.results_tree.heading('Score', text='Score')
        self.results_tree.heading('Warnings', text='Warnings')
        self.results_tree.heading('Errors', text='Errors')
        
        # Column widths
        self.results_tree.column('Rank', width=50, anchor=tk.CENTER)
        self.results_tree.column('Vector Fwd', width=100, anchor=tk.CENTER)
        self.results_tree.column('Vector Rev', width=100, anchor=tk.CENTER)
        self.results_tree.column('Insert Fwd', width=100, anchor=tk.CENTER)
        self.results_tree.column('Insert Rev', width=100, anchor=tk.CENTER)
        self.results_tree.column('Score', width=80, anchor=tk.CENTER)
        self.results_tree.column('Warnings', width=80, anchor=tk.CENTER)
        self.results_tree.column('Errors', width=80, anchor=tk.CENTER)
        
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
        ttk.Button(button_frame, text="Export All", 
                  command=self.export_all).pack(side=tk.LEFT, padx=5)
    
    def setup_details_tab(self, parent):
        """Setup primer details tab"""
        
        # Details text
        self.details_text = scrolledtext.ScrolledText(parent, height=35, width=100, 
                                                      wrap=tk.WORD, font=('Courier', 9))
        self.details_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure tags for formatting
        self.details_text.tag_configure('header', font=('Courier', 10, 'bold'))
        self.details_text.tag_configure('subheader', font=('Courier', 9, 'bold'))
        self.details_text.tag_configure('good', foreground='green')
        self.details_text.tag_configure('warning', foreground='orange')
        self.details_text.tag_configure('error', foreground='red')
        self.details_text.tag_configure('highlight', foreground='purple', font=('Courier', 9, 'bold'))
    
    def toggle_insert_range(self):
        """Enable/disable insert range inputs"""
        if self.optimize_insert_var.get():
            self.insert_fwd_start_entry.config(state='normal')
            self.insert_fwd_end_entry.config(state='normal')
            self.insert_rev_start_entry.config(state='normal')
            self.insert_rev_end_entry.config(state='normal')
        else:
            self.insert_fwd_start_entry.config(state='disabled')
            self.insert_fwd_end_entry.config(state='disabled')
            self.insert_rev_start_entry.config(state='disabled')
            self.insert_rev_end_entry.config(state='disabled')
    
    def update_vector_length(self, event=None):
        """Update vector sequence length display"""
        seq = self.vector_text.get(1.0, tk.END).strip().replace('\n', '').replace(' ', '')
        self.vector_length_var.set(f"Length: {len(seq)} bp")
    
    def update_insert_length(self, event=None):
        """Update insert sequence length display"""
        seq = self.insert_text.get(1.0, tk.END).strip().replace('\n', '').replace(' ', '')
        self.insert_length_var.set(f"Length: {len(seq)} bp")
    
    def load_vector_file(self):
        """Load vector sequence from file"""
        filename = filedialog.askopenfilename(
            title="Select Vector File",
            filetypes=[("All Files", "*.*"), ("FASTA", "*.fasta *.fa"), 
                      ("GenBank", "*.gb *.gbk"), ("Text", "*.txt")]
        )
        if filename:
            try:
                from utils import read_sequence_file
                sequence = read_sequence_file(filename)
                self.vector_text.delete(1.0, tk.END)
                self.vector_text.insert(1.0, sequence)
                # Store for GenBank export
                self.vector_file = filename if filename.lower().endswith(('.gb', '.gbk', '.genbank')) else None
                self.vector_seq = sequence
                # Update length display
                self.update_vector_length()
                self.status_var.set(f"Loaded vector: {len(sequence)} bp")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {e}")
    
    def load_insert_file(self):
        """Load insert sequence from file"""
        filename = filedialog.askopenfilename(
            title="Select Insert File",
            filetypes=[("All Files", "*.*"), ("FASTA", "*.fasta *.fa"), 
                      ("GenBank", "*.gb *.gbk"), ("Text", "*.txt")]
        )
        if filename:
            try:
                from utils import read_sequence_file
                sequence = read_sequence_file(filename)
                self.insert_text.delete(1.0, tk.END)
                self.insert_text.insert(1.0, sequence)
                # Store for GenBank export
                self.insert_file = filename if filename.lower().endswith(('.gb', '.gbk', '.genbank')) else None
                self.insert_seq = sequence
                # Update length display
                self.update_insert_length()
                self.status_var.set(f"Loaded insert: {len(sequence)} bp")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to load file: {e}")
    
    def clear_all(self):
        """Clear all inputs"""
        self.vector_text.delete(1.0, tk.END)
        self.insert_text.delete(1.0, tk.END)
        self.results_tree.delete(*self.results_tree.get_children())
        self.details_text.delete(1.0, tk.END)
        self.current_results = []
        self.vector_file = None
        self.insert_file = None
        self.vector_seq = None
        self.insert_seq = None
        self.status_var.set("Cleared all inputs")
    
    def run_optimization(self):
        """Run primer optimization in background thread"""
        
        # Get sequences
        vector_seq = self.vector_text.get(1.0, tk.END).strip().replace('\n', '').replace(' ', '').upper()
        insert_seq = self.insert_text.get(1.0, tk.END).strip().replace('\n', '').replace(' ', '').upper()
        
        # Store sequences for GenBank export (if not already stored from file load)
        if not self.vector_seq:
            self.vector_seq = vector_seq
        if not self.insert_seq:
            self.insert_seq = insert_seq
        
        # Validate inputs
        if not vector_seq or not insert_seq:
            messagebox.showerror("Error", "Please enter both vector and insert sequences")
            return
        
        if len(vector_seq) < 100:
            messagebox.showerror("Error", "Vector sequence too short (minimum 100 bp)")
            return
        
        if len(insert_seq) < 20:
            messagebox.showerror("Error", "Insert sequence too short (minimum 20 bp)")
            return
        
        # Get parameters
        try:
            vector_upstream_start = int(self.vector_upstream_start_var.get())
            vector_upstream_end = int(self.vector_upstream_end_var.get())
            vector_downstream_start = int(self.vector_downstream_start_var.get())
            vector_downstream_end = int(self.vector_downstream_end_var.get())
            homology_length = int(self.homology_var.get())
            n_results = int(self.n_results_var.get())
            step = int(self.step_var.get())
            
            # Get insert range parameters if optimization is enabled
            optimize_insert = self.optimize_insert_var.get()
            if optimize_insert:
                insert_fwd_start = int(self.insert_fwd_start_var.get())
                insert_fwd_end = int(self.insert_fwd_end_var.get())
                insert_rev_start = int(self.insert_rev_start_var.get())
                insert_rev_end = int(self.insert_rev_end_var.get())
            else:
                insert_fwd_start = insert_fwd_end = insert_rev_start = insert_rev_end = None
        except ValueError:
            messagebox.showerror("Error", "Invalid parameter values")
            return
        
        # Validate vector ranges
        if vector_upstream_start >= vector_upstream_end:
            messagebox.showerror("Error", "Vector forward primer start must be less than end")
            return
        
        if vector_downstream_start >= vector_downstream_end:
            messagebox.showerror("Error", "Vector reverse primer start must be less than end")
            return
        
        if vector_upstream_end >= vector_downstream_start:
            messagebox.showerror("Error", "Vector forward primer region must not overlap with reverse primer region")
            return
        
        if vector_downstream_end > len(vector_seq):
            messagebox.showerror("Error", f"Vector reverse primer end exceeds vector length ({len(vector_seq)} bp)")
            return
        
        # Validate insert ranges if optimization is enabled
        if optimize_insert:
            if insert_fwd_start >= insert_fwd_end:
                messagebox.showerror("Error", "Insert forward primer start must be less than end")
                return
            
            if insert_rev_start >= insert_rev_end:
                messagebox.showerror("Error", "Insert reverse primer start must be less than end")
                return
            
            if insert_fwd_end >= insert_rev_start:
                messagebox.showerror("Error", "Insert forward primer region must not overlap with reverse primer region")
                return
            
            if insert_rev_end > len(insert_seq):
                messagebox.showerror("Error", f"Insert reverse primer end exceeds insert length ({len(insert_seq)} bp)")
                return
        
        # Warn about complexity
        vec_up_tests = ((vector_upstream_end - vector_upstream_start) // step + 1)
        vec_down_tests = ((vector_downstream_end - vector_downstream_start) // step + 1)
        total_tests = vec_up_tests * vec_down_tests
        
        if optimize_insert:
            ins_fwd_tests = ((insert_fwd_end - insert_fwd_start) // step + 1)
            ins_rev_tests = ((insert_rev_end - insert_rev_start) // step + 1)
            total_tests *= ins_fwd_tests * ins_rev_tests
        
        if total_tests > 50000:
            result = messagebox.askyesno(
                "High Complexity Warning",
                f"This will test {total_tests:,} combinations and may take 5-10 minutes!\n\n"
                f"Recommendations:\n"
                f"• Increase step size (currently {step} bp)\n"
                f"• Reduce range sizes\n"
                f"• Disable insert optimization\n\n"
                f"Continue anyway?"
            )
            if not result:
                return
        elif total_tests > 10000:
            messagebox.showinfo(
                "Complexity Notice",
                f"Testing {total_tests:,} combinations. This may take 1-2 minutes.\n\n"
                f"Tip: Increase step size for faster results."
            )
        
        # Disable button and start progress
        self.optimize_btn.config(state='disabled')
        self.progress.start()
        self.status_var.set("Optimizing primer design...")
        
        # Run in thread
        def optimize_thread():
            try:
                results = self.optimizer.get_top_n_results(
                    vector_seq, insert_seq, 
                    vector_upstream_start=vector_upstream_start,
                    vector_upstream_end=vector_upstream_end,
                    vector_downstream_start=vector_downstream_start,
                    vector_downstream_end=vector_downstream_end,
                    insert_fwd_start=insert_fwd_start,
                    insert_fwd_end=insert_fwd_end,
                    insert_rev_start=insert_rev_start,
                    insert_rev_end=insert_rev_end,
                    optimize_insert=optimize_insert,
                    n=n_results, 
                    homology_length=homology_length,
                    step=step
                )
                
                # Update UI in main thread
                self.root.after(0, lambda: self.display_results(results))
                
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror("Error", f"Optimization failed: {e}"))
                self.root.after(0, lambda: self.status_var.set("Optimization failed"))
            finally:
                self.root.after(0, lambda: self.progress.stop())
                self.root.after(0, lambda: self.optimize_btn.config(state='normal'))
        
        thread = threading.Thread(target=optimize_thread, daemon=True)
        thread.start()
    
    def display_results(self, results):
        """Display optimization results"""
        
        # Clear existing results
        self.results_tree.delete(*self.results_tree.get_children())
        self.current_results = results
        
        # Populate tree
        for i, result in enumerate(results, 1):
            self.results_tree.insert('', tk.END, values=(
                i,
                result.get('vector_upstream_pos', '-'),
                result.get('vector_downstream_pos', '-'),
                result.get('insert_fwd_pos', '-'),
                result.get('insert_rev_pos', '-'),
                f"{result['score']:.1f}",
                result['warnings'],
                result['errors']
            ))
        
        self.status_var.set(f"Found {len(results)} optimal combinations")
        
        # Select first result
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
    
    def view_selected_details(self):
        """View details of selected result"""
        selection = self.results_tree.selection()
        if not selection:
            messagebox.showinfo("Info", "Please select a result to view details")
            return
        
        # Switch to details tab
        # (Assuming notebook is accessible - you might need to store reference)
    
    def display_primer_details(self, result):
        """Display detailed primer information"""
        
        self.details_text.delete(1.0, tk.END)
        
        # Header
        self.details_text.insert(tk.END, "="*80 + "\n", 'header')
        self.details_text.insert(tk.END, "PRIMER DESIGN DETAILS\n", 'header')
        self.details_text.insert(tk.END, f"Vector Forward Primer Position: {result.get('vector_upstream_pos', '-')} bp\n", 'header')
        self.details_text.insert(tk.END, f"Vector Reverse Primer Position: {result.get('vector_downstream_pos', '-')} bp\n", 'header')
        self.details_text.insert(tk.END, f"Insert Forward Primer Position: {result.get('insert_fwd_pos', '-')} bp\n", 'header')
        self.details_text.insert(tk.END, f"Insert Reverse Primer Position: {result.get('insert_rev_pos', '-')} bp\n", 'header')
        self.details_text.insert(tk.END, f"Overall Score: {result['score']:.1f}/100\n", 'header')
        self.details_text.insert(tk.END, "="*80 + "\n\n", 'header')
        
        primers = result['primers']
        vals = result['validations']
        
        # Vector primers
        self.details_text.insert(tk.END, "VECTOR LINEARIZATION PRIMERS\n", 'subheader')
        self.details_text.insert(tk.END, "-"*80 + "\n")
        
        self.details_text.insert(tk.END, "\n1. Vector Forward:\n", 'subheader')
        self.display_formatted_primer(primers['vector_forward'], vals['vector_forward'])
        self.display_validation_summary(vals['vector_forward'])
        
        self.details_text.insert(tk.END, "\n2. Vector Reverse:\n", 'subheader')
        self.display_formatted_primer(primers['vector_reverse'], vals['vector_reverse'])
        self.display_validation_summary(vals['vector_reverse'])
        
        # Insert primers
        self.details_text.insert(tk.END, "\n\nINSERT AMPLIFICATION PRIMERS\n", 'subheader')
        self.details_text.insert(tk.END, "-"*80 + "\n")
        
        self.details_text.insert(tk.END, "\n3. Insert Forward:\n", 'subheader')
        self.display_formatted_primer(primers['insert_forward'], vals['insert_forward'])
        self.display_validation_summary(vals['insert_forward'])
        
        self.details_text.insert(tk.END, "\n4. Insert Reverse:\n", 'subheader')
        self.display_formatted_primer(primers['insert_reverse'], vals['insert_reverse'])
        self.display_validation_summary(vals['insert_reverse'])
        
        # Final construct
        if 'final_construct' in result:
            self.details_text.insert(tk.END, "\n\nFINAL ASSEMBLED CONSTRUCT\n", 'subheader')
            self.details_text.insert(tk.END, "-"*80 + "\n")
            final_seq = result['final_construct']
            self.details_text.insert(tk.END, f"Length: {len(final_seq)} bp\n", 'good')
            self.details_text.insert(tk.END, f"\nSequence (5' to 3'):\n")
            self.details_text.insert(tk.END, "Note: Homology regions shown in UPPERCASE\n\n", 'good')
            
            # Get homology arm lengths from validations
            homology_len_fwd = vals.get('insert_forward', {}).get('homology_region', {}).get('length', 25)
            homology_len_rev = vals.get('insert_reverse', {}).get('homology_region', {}).get('length', 25)
            
            # Calculate positions for highlighting
            vector_upstream_len = result.get('vector_upstream_pos', 0)
            insert_len = result.get('insert_rev_pos', 0) - result.get('insert_fwd_pos', 0)
            
            # Homology regions at junctions
            homology1_start = vector_upstream_len - homology_len_fwd
            homology1_end = vector_upstream_len + homology_len_fwd
            homology2_start = vector_upstream_len + insert_len - homology_len_rev
            homology2_end = vector_upstream_len + insert_len + homology_len_rev
            
            # Display sequence in chunks of 60 bp with highlighting
            for i in range(0, len(final_seq), 60):
                chunk_start = i
                chunk_end = min(i + 60, len(final_seq))
                line_num = f"{i+1:6d}  "
                self.details_text.insert(tk.END, line_num)
                
                # Process each character in the chunk
                for pos in range(chunk_start, chunk_end):
                    char = final_seq[pos]
                    # Check if position is in homology region
                    if (homology1_start <= pos < homology1_end) or (homology2_start <= pos < homology2_end):
                        self.details_text.insert(tk.END, char.upper(), 'highlight')
                    else:
                        self.details_text.insert(tk.END, char.lower())
                
                self.details_text.insert(tk.END, "\n")
    
    def display_formatted_primer(self, primer_seq, validation):
        """Display primer with homology in lowercase and annealing in uppercase"""
        if 'homology_region' in validation and 'annealing_region' in validation:
            homology_seq = validation['homology_region']['sequence']
            annealing_seq = validation['annealing_region']['sequence']
            # Format: homology (lowercase) + annealing (uppercase)
            formatted = homology_seq.lower() + annealing_seq.upper()
            self.details_text.insert(tk.END, f"   5'- ")
            self.details_text.insert(tk.END, homology_seq.lower())
            self.details_text.insert(tk.END, annealing_seq.upper(), 'highlight')
            self.details_text.insert(tk.END, " -3'\n")
        else:
            self.details_text.insert(tk.END, f"   5'- {primer_seq} -3'\n")
    
    def display_validation_summary(self, validation):
        """Display validation summary for a primer"""
        
        metrics = validation['metrics']
        
        self.details_text.insert(tk.END, f"   Length: {metrics['length']} nt\n")
        self.details_text.insert(tk.END, f"   GC: {metrics['gc_content']}%\n")
        
        if 'annealing_region' in validation:
            self.details_text.insert(tk.END, f"   Tm (annealing): {validation['annealing_region']['tm']}°C\n")
            if 'homology_region' in validation:
                self.details_text.insert(tk.END, f"   Tm (homology): {validation['homology_region']['tm']:.1f}°C\n")
        else:
            self.details_text.insert(tk.END, f"   Tm: {metrics['tm']}°C\n")
        
        # Warnings
        if validation['warnings']:
            self.details_text.insert(tk.END, f"   ⚠ Warnings: {len(validation['warnings'])}\n", 'warning')
            for warning in validation['warnings']:
                self.details_text.insert(tk.END, f"      • {warning}\n", 'warning')
        else:
            self.details_text.insert(tk.END, "   ✓ No warnings\n", 'good')
        
        # Errors
        if validation['errors']:
            self.details_text.insert(tk.END, f"   ✗ Errors: {len(validation['errors'])}\n", 'error')
            for error in validation['errors']:
                self.details_text.insert(tk.END, f"      • {error}\n", 'error')
    
    def export_selected(self):
        """Export selected result to file"""
        selection = self.results_tree.selection()
        if not selection:
            messagebox.showinfo("Info", "Please select a result to export")
            return
        
        item = self.results_tree.item(selection[0])
        rank = int(item['values'][0]) - 1
        
        if 0 <= rank < len(self.current_results):
            self.show_export_dialog(self.current_results[rank])
    
    def export_all(self):
        """Export all results to file"""
        if not self.current_results:
            messagebox.showinfo("Info", "No results to export")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Export All Results",
            defaultextension=".txt",
            filetypes=[("Text", "*.txt"), ("All Files", "*.*")]
        )
        
        if filename:
            try:
                with open(filename, 'w') as f:
                    f.write("Gibson Assembly Primer Design - All Results\n")
                    f.write("="*80 + "\n\n")
                    
                    for i, result in enumerate(self.current_results, 1):
                        f.write(f"\n{'='*80}\n")
                        f.write(f"RESULT #{i} - Cut Site: {result['cut_site']} bp\n")
                        f.write(f"Score: {result['score']:.1f}/100\n")
                        f.write(f"{'='*80}\n\n")
                        
                        # Write primers
                        f.write("PRIMERS:\n")
                        f.write(f"Vector Forward:  5'- {result['primers']['vector_forward']} -3'\n")
                        f.write(f"Vector Reverse:  5'- {result['primers']['vector_reverse']} -3'\n")
                        f.write(f"Insert Forward:  5'- {result['primers']['insert_forward']} -3'\n")
                        f.write(f"Insert Reverse:  5'- {result['primers']['insert_reverse']} -3'\n\n")
                        
                        # Write final construct
                        if 'final_construct' in result:
                            f.write("FINAL CONSTRUCT:\n")
                            f.write(f"Length: {len(result['final_construct'])} bp\n")
                            f.write(f"Sequence: {result['final_construct']}\n\n")
                
                messagebox.showinfo("Success", f"Exported all results to {filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {e}")
    
    def show_export_dialog(self, result):
        """Show export options dialog"""
        dialog = tk.Toplevel(self.root)
        dialog.title("Export Options")
        dialog.geometry("450x350")
        dialog.transient(self.root)
        dialog.grab_set()
        
        # Center the dialog
        dialog.update_idletasks()
        x = (dialog.winfo_screenwidth() // 2) - (450 // 2)
        y = (dialog.winfo_screenheight() // 2) - (350 // 2)
        dialog.geometry(f"450x350+{x}+{y}")
        
        main_frame = ttk.Frame(dialog, padding="20")
        main_frame.pack(fill=tk.BOTH, expand=True)
        
        # Title
        ttk.Label(main_frame, text="Choose Export Options", font=('Arial', 12, 'bold')).pack(pady=(0, 15))
        
        # Export content selection
        content_frame = ttk.LabelFrame(main_frame, text="What to Export", padding="10")
        content_frame.pack(fill=tk.X, pady=(0, 10))
        
        export_content = tk.StringVar(value="construct")
        ttk.Radiobutton(content_frame, text="Final Construct Only", variable=export_content, value="construct").pack(anchor=tk.W, pady=2)
        ttk.Radiobutton(content_frame, text="Primers Only", variable=export_content, value="primers").pack(anchor=tk.W, pady=2)
        ttk.Radiobutton(content_frame, text="Both Construct and Primers", variable=export_content, value="both").pack(anchor=tk.W, pady=2)
        
        # Format selection
        format_frame = ttk.LabelFrame(main_frame, text="Export Format", padding="10")
        format_frame.pack(fill=tk.X, pady=(0, 15))
        
        export_format = tk.StringVar(value="genbank")
        ttk.Radiobutton(format_frame, text="GenBank (.gb) - with annotations", variable=export_format, value="genbank").pack(anchor=tk.W, pady=2)
        ttk.Radiobutton(format_frame, text="FASTA (.fasta) - sequence only", variable=export_format, value="fasta").pack(anchor=tk.W, pady=2)
        ttk.Radiobutton(format_frame, text="Text (.txt) - detailed report", variable=export_format, value="text").pack(anchor=tk.W, pady=2)
        
        # Buttons
        button_frame = ttk.Frame(main_frame)
        button_frame.pack(pady=(10, 0))
        
        def on_export():
            dialog.destroy()
            self.export_result(result, export_content.get(), export_format.get())
        
        ttk.Button(button_frame, text="Export", command=on_export, width=12).pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Cancel", command=dialog.destroy, width=12).pack(side=tk.LEFT, padx=5)
    
    def export_result(self, result, content_type="both", format_type="genbank"):
        """Export single result to file"""
        # Determine file extension and filter
        if format_type == "genbank":
            default_ext = ".gb"
            filetypes = [("GenBank", "*.gb *.gbk"), ("All Files", "*.*")]
        elif format_type == "fasta":
            default_ext = ".fasta"
            filetypes = [("FASTA", "*.fasta *.fa"), ("All Files", "*.*")]
        else:
            default_ext = ".txt"
            filetypes = [("Text", "*.txt"), ("All Files", "*.*")]
        
        filename = filedialog.asksaveasfilename(
            title="Export Result",
            defaultextension=default_ext,
            filetypes=filetypes
        )
        
        if filename:
            try:
                if format_type == "genbank":
                    # Export as GenBank with annotations
                    from utils import create_genbank_output, save_genbank
                    from Bio.SeqRecord import SeqRecord
                    from Bio.Seq import Seq
                    from Bio.SeqFeature import SeqFeature, FeatureLocation
                    
                    if content_type == "construct":
                        # Export only construct
                        record = create_genbank_output(
                            result,
                            vector_seq=self.vector_seq,
                            insert_seq=self.insert_seq,
                            original_vector_file=self.vector_file,
                            original_insert_file=self.insert_file
                        )
                        save_genbank(record, filename)
                    
                    elif content_type == "primers":
                        # Export primers as separate sequences
                        from Bio import SeqIO
                        primers = result['primers']
                        records = []
                        for name, seq in primers.items():
                            rec = SeqRecord(
                                Seq(seq),
                                id=name.replace('_', '_'),
                                name=name,
                                description=f"Gibson assembly primer - {name}"
                            )
                            records.append(rec)
                        with open(filename, 'w') as f:
                            SeqIO.write(records, f, "genbank")
                    
                    else:  # both
                        record = create_genbank_output(
                            result,
                            vector_seq=self.vector_seq,
                            insert_seq=self.insert_seq,
                            original_vector_file=self.vector_file,
                            original_insert_file=self.insert_file
                        )
                        save_genbank(record, filename)
                    
                    messagebox.showinfo("Success", f"Exported GenBank file to:\n{filename}")
                
                elif format_type == "fasta":
                    # Export as FASTA
                    from Bio.SeqRecord import SeqRecord
                    from Bio.Seq import Seq
                    from Bio import SeqIO
                    
                    records = []
                    
                    if content_type in ["construct", "both"]:
                        record = SeqRecord(
                            Seq(result['final_construct']),
                            id="gibson_construct",
                            description=f"Gibson assembly - Score: {result['score']:.1f}/100"
                        )
                        records.append(record)
                    
                    if content_type in ["primers", "both"]:
                        primers = result['primers']
                        for name, seq in primers.items():
                            rec = SeqRecord(
                                Seq(seq),
                                id=name,
                                description=f"Gibson assembly primer"
                            )
                            records.append(rec)
                    
                    with open(filename, 'w') as f:
                        SeqIO.write(records, f, "fasta")
                    messagebox.showinfo("Success", f"Exported FASTA to:\n{filename}")
                
                else:  # text format
                    # Export as plain text
                    with open(filename, 'w') as f:
                        f.write("Gibson Assembly Primer Design\n")
                        f.write("="*80 + "\n\n")
                        f.write(f"Vector Forward Primer Position: {result.get('vector_upstream_pos', '-')} bp\n")
                        f.write(f"Vector Reverse Primer Position: {result.get('vector_downstream_pos', '-')} bp\n")
                        f.write(f"Insert Forward Primer Position: {result.get('insert_fwd_pos', '-')} bp\n")
                        f.write(f"Insert Reverse Primer Position: {result.get('insert_rev_pos', '-')} bp\n")
                        f.write(f"Overall Score: {result['score']:.1f}/100\n\n")
                        
                        if content_type in ["primers", "both"]:
                            f.write("VECTOR LINEARIZATION PRIMERS\n")
                            f.write("-"*80 + "\n")
                            for primer_name, label in [('vector_forward', 'Forward'), ('vector_reverse', 'Reverse')]:
                                f.write(f"\n{label}: 5'- {result['primers'][primer_name]} -3'\n")
                                val = result['validations'][primer_name]
                                f.write(f"  Length: {val['metrics']['length']} nt\n")
                                f.write(f"  GC: {val['metrics']['gc_content']}%\n")
                                if 'annealing_region' in val:
                                    f.write(f"  Tm (annealing): {val['annealing_region']['tm']}°C\n")
                                    if 'homology_region' in val:
                                        f.write(f"  Tm (homology): {val['homology_region']['tm']:.1f}°C\n")
                                else:
                                    f.write(f"  Tm: {val['metrics']['tm']}°C\n")
                            
                            f.write("\n\nINSERT AMPLIFICATION PRIMERS\n")
                            f.write("-"*80 + "\n")
                            for primer_name, label in [('insert_forward', 'Forward'), ('insert_reverse', 'Reverse')]:
                                f.write(f"\n{label}: 5'- {result['primers'][primer_name]} -3'\n")
                                val = result['validations'][primer_name]
                                f.write(f"  Length: {val['metrics']['length']} nt\n")
                                f.write(f"  GC: {val['metrics']['gc_content']}%\n")
                                if 'annealing_region' in val:
                                    f.write(f"  Tm (annealing): {val['annealing_region']['tm']}°C\n")
                                    if 'homology_region' in val:
                                        f.write(f"  Tm (homology): {val['homology_region']['tm']:.1f}°C\n")
                                else:
                                    f.write(f"  Tm: {val['metrics']['tm']}°C\n")
                            f.write("\n")
                        
                        # Write final construct
                        if content_type in ["construct", "both"] and 'final_construct' in result:
                            f.write("\nFINAL ASSEMBLED CONSTRUCT\n")
                            f.write("-"*80 + "\n")
                            f.write(f"Length: {len(result['final_construct'])} bp\n\n")
                            f.write("Sequence (5' to 3'):\n")
                            final_seq = result['final_construct']
                            for i in range(0, len(final_seq), 60):
                                chunk = final_seq[i:i+60]
                                f.write(f"{i+1:6d}  {chunk}\n")
                    
                    messagebox.showinfo("Success", f"Exported to:\n{filename}")
            
            except Exception as e:
                messagebox.showerror("Error", f"Export failed: {e}")


def main():
    """Main entry point for GUI"""
    root = tk.Tk()
    app = GibsonAssemblyGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
