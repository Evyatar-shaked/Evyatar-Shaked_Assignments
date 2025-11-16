"""
BLAST Search GUI Application
A graphical interface for performing NCBI BLAST searches on DNA sequences
"""

import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import threading
from blast_search import BLASTSearcher


class BLASTApp:
    """Main GUI application for BLAST searches"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("DNA BLAST Search Tool")
        self.root.geometry("900x700")
        
        self.searcher = BLASTSearcher()
        self.search_thread = None
        
        self.setup_ui()
        
    def setup_ui(self):
        """Setup the user interface"""
        
        # Main frame with padding
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(2, weight=1)
        
        # Title
        title_label = ttk.Label(
            main_frame, 
            text="NCBI BLAST DNA Sequence Search", 
            font=("Arial", 16, "bold")
        )
        title_label.grid(row=0, column=0, pady=10)
        
        # Input frame
        input_frame = ttk.LabelFrame(main_frame, text="Input Parameters", padding="10")
        input_frame.grid(row=1, column=0, sticky=(tk.W, tk.E), pady=5)
        input_frame.columnconfigure(1, weight=1)
        
        # DNA Sequence input
        ttk.Label(input_frame, text="DNA Sequence:").grid(row=0, column=0, sticky=tk.W, pady=5)
        
        self.sequence_text = scrolledtext.ScrolledText(
            input_frame, 
            height=8, 
            width=70,
            wrap=tk.WORD,
            font=("Courier", 10)
        )
        self.sequence_text.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        
        # Organism selection
        organism_frame = ttk.Frame(input_frame)
        organism_frame.grid(row=2, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=10)
        organism_frame.columnconfigure(1, weight=1)
        
        ttk.Label(organism_frame, text="Organism:").grid(row=0, column=0, sticky=tk.W, padx=(0, 10))
        
        # Organism dropdown
        self.organism_var = tk.StringVar()
        self.organism_combo = ttk.Combobox(
            organism_frame, 
            textvariable=self.organism_var,
            width=40,
            state='normal'
        )
        
        # Common organisms list
        common_organisms = [
            "All organisms",
            "Homo sapiens",
            "Mus musculus",
            "Rattus norvegicus",
            "Drosophila melanogaster",
            "Caenorhabditis elegans",
            "Saccharomyces cerevisiae",
            "Escherichia coli",
            "Arabidopsis thaliana",
            "Danio rerio",
            "Xenopus laevis"
        ]
        
        self.organism_combo['values'] = common_organisms
        self.organism_combo.current(0)
        self.organism_combo.grid(row=0, column=1, sticky=(tk.W, tk.E))
        
        # Max results
        ttk.Label(organism_frame, text="Max Results:").grid(row=0, column=2, sticky=tk.W, padx=(20, 10))
        
        self.max_results_var = tk.StringVar(value="10")
        max_results_spinbox = ttk.Spinbox(
            organism_frame,
            from_=1,
            to=50,
            textvariable=self.max_results_var,
            width=10
        )
        max_results_spinbox.grid(row=0, column=3, sticky=tk.W)
        
        # Search button
        self.search_button = ttk.Button(
            input_frame,
            text="Search BLAST",
            command=self.perform_search,
            style='Accent.TButton'
        )
        self.search_button.grid(row=3, column=0, columnspan=2, pady=10)
        
        # Progress bar
        self.progress = ttk.Progressbar(
            input_frame,
            mode='indeterminate',
            length=400
        )
        self.progress.grid(row=4, column=0, columnspan=2, pady=5)
        
        # Status label
        self.status_label = ttk.Label(
            input_frame,
            text="Ready to search",
            foreground="green"
        )
        self.status_label.grid(row=5, column=0, columnspan=2)
        
        # Results frame
        results_frame = ttk.LabelFrame(main_frame, text="Search Results", padding="10")
        results_frame.grid(row=2, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)
        
        # Results text area
        self.results_text = scrolledtext.ScrolledText(
            results_frame,
            height=20,
            width=80,
            wrap=tk.WORD,
            font=("Courier", 9)
        )
        self.results_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Clear button
        clear_button = ttk.Button(
            results_frame,
            text="Clear Results",
            command=self.clear_results
        )
        clear_button.grid(row=1, column=0, pady=5)
        
        # Example sequence button
        example_button = ttk.Button(
            input_frame,
            text="Load Example Sequence",
            command=self.load_example_sequence
        )
        example_button.grid(row=6, column=0, columnspan=2, pady=5)
        
    def load_example_sequence(self):
        """Load an example DNA sequence"""
        example_seq = """ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCCGCCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC"""
        
        self.sequence_text.delete('1.0', tk.END)
        self.sequence_text.insert('1.0', example_seq)
        self.organism_var.set("Homo sapiens")
        
        messagebox.showinfo(
            "Example Loaded",
            "Loaded a fragment of the human insulin gene (INS)"
        )
    
    def perform_search(self):
        """Perform BLAST search in a separate thread"""
        
        # Get sequence
        sequence = self.sequence_text.get('1.0', tk.END).strip()
        
        if not sequence:
            messagebox.showerror("Error", "Please enter a DNA sequence")
            return
        
        # Validate sequence
        is_valid, result = self.searcher.validate_sequence(sequence)
        
        if not is_valid:
            messagebox.showerror("Invalid Sequence", result)
            return
        
        # Get organism
        organism = self.organism_var.get()
        if organism == "All organisms":
            organism = None
        
        # Get max results
        try:
            max_results = int(self.max_results_var.get())
        except ValueError:
            max_results = 10
        
        # Disable search button and start progress
        self.search_button.config(state='disabled')
        self.progress.start(10)
        self.status_label.config(text="Searching NCBI BLAST database...", foreground="blue")
        self.results_text.delete('1.0', tk.END)
        self.results_text.insert('1.0', "Submitting query to NCBI BLAST...\nThis may take 1-2 minutes...\n\n")
        
        # Run search in thread
        self.search_thread = threading.Thread(
            target=self._search_worker,
            args=(result, organism, max_results),
            daemon=True
        )
        self.search_thread.start()
    
    def _search_worker(self, sequence, organism, max_results):
        """Worker thread for BLAST search"""
        try:
            results = self.searcher.search_sequence(sequence, organism, max_results)
            self.root.after(0, self._display_results, results)
        except Exception as e:
            self.root.after(0, self._display_error, str(e))
    
    def _display_results(self, results):
        """Display search results in the text area"""
        self.progress.stop()
        self.search_button.config(state='normal')
        
        if not results:
            self.status_label.config(text="No results found", foreground="orange")
            self.results_text.delete('1.0', tk.END)
            self.results_text.insert('1.0', "No matching sequences found in the database.")
            return
        
        self.status_label.config(text=f"Search complete - {len(results)} results found", foreground="green")
        
        # Display results
        self.results_text.delete('1.0', tk.END)
        
        header = f"{'='*80}\n"
        header += f"BLAST SEARCH RESULTS - {len(results)} matches found\n"
        header += f"{'='*80}\n\n"
        
        self.results_text.insert('1.0', header)
        
        for i, result in enumerate(results, 1):
            result_text = f"\n{'='*80}\n"
            result_text += f"RESULT #{i}\n"
            result_text += f"{'='*80}\n"
            result_text += f"Title: {result['title']}\n"
            result_text += f"Accession: {result['accession']}\n"
            result_text += f"Sequence Length: {result['length']} bp\n"
            result_text += f"\n--- Alignment Statistics ---\n"
            result_text += f"Identity: {result['identity_percent']:.2f}%\n"
            result_text += f"E-value: {result['e_value']:.2e}\n"
            result_text += f"Score: {result['score']}\n"
            result_text += f"Identities: {result['identities']}/{result['alignment_length']}\n"
            result_text += f"Gaps: {result['gaps']}/{result['alignment_length']}\n"
            result_text += f"\n--- Alignment Preview (first 200 characters) ---\n"
            result_text += f"Query:   {result['query'][:200]}\n"
            result_text += f"Match:   {result['match'][:200]}\n"
            result_text += f"Subject: {result['subject'][:200]}\n"
            result_text += "\n"
            
            self.results_text.insert(tk.END, result_text)
        
        # Scroll to top
        self.results_text.see('1.0')
    
    def _display_error(self, error_msg):
        """Display error message"""
        self.progress.stop()
        self.search_button.config(state='normal')
        self.status_label.config(text="Search failed", foreground="red")
        
        messagebox.showerror("Search Error", f"BLAST search failed:\n\n{error_msg}")
    
    def clear_results(self):
        """Clear the results text area"""
        self.results_text.delete('1.0', tk.END)
        self.status_label.config(text="Ready to search", foreground="green")


def main():
    """Main entry point for the application"""
    root = tk.Tk()
    app = BLASTApp(root)
    root.mainloop()


if __name__ == "__main__":
    main()
