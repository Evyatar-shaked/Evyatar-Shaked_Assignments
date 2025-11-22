# Gibson Assembly Primer Design Wizard 🧬

A Python-based tool for designing primers for **binary Gibson assembly** following the **Takara InFusion protocol** (50°C incubation). This wizard automates the design of primers for vector linearization and insert amplification with optimal homology regions.

## Features

- ✅ **Binary Assembly Support**: Designs primers for vector linearization + insert amplification
- ✅ **Graphical Interface (GUI)**: Modern tkinter-based interface for easy use
- ✅ **Automatic Cut Site Optimization**: Scans vector range to find best primer designs
- ✅ **Biopython Integration**: Uses industry-standard Biopython library for accurate calculations
- ✅ **Intelligent Scoring System**: Ranks primers by Tm, GC content, dimers, and more
- ✅ **Accurate Tm Calculation**: Nearest-neighbor thermodynamics with proper salt corrections
- ✅ **Comprehensive Validation**: Checks for hairpins, self-dimers, hetero-dimers, and secondary structures
- ✅ **Takara InFusion Optimized**: Follows recommended parameters for InFusion assembly
- ✅ **Interactive CLI Wizard**: Command-line interface for terminal users
- ✅ **Multiple Results**: View top N best cut sites with detailed comparisons
- ✅ **Detailed Reports**: Generates validation reports and primer order sheets
- ✅ **Multiple Format Support**: Reads FASTA, GenBank, and plain text files

## Installation

1. Clone or download this repository
2. Install Biopython:

```bash
pip install -r requirements.txt
```

Or manually:
```bash
pip install biopython
```

3. Run the wizard:
```bash
cd gibson_assembly2.0
python gibson_wizard.py
```

## Quick Start

**GUI Application (Recommended):**

```bash
python gibson_gui.py
```

**Command-line Wizard:**

```bash
python gibson_wizard.py
```

The wizard will guide you through:
1. Entering your vector (plasmid) sequence
2. Specifying the linearization/cut site
3. Entering your insert (gene) sequence
4. Designing homology regions
5. Generating optimized primers with validation

## Design Rules

### Standard PCR Primer Rules (Annealing Region)

- **Length**: 18-25 nt
- **Tm**: 55-65°C (calculated using nearest-neighbor method)
- **GC Content**: 40-60%
- **GC Clamp**: 1-2 G/C bases at 3' end preferred
- **Tm Difference**: ≤3°C between forward and reverse primers

### Gibson-Specific Rules (Homology Overhang)

- **Homology Length**: 15-40 nt (recommended: 20-30 nt)
- **Homology Tm**: ~50-60°C for stability
- **Sequence Design**: 
  - Unique to adjacent fragment
  - Avoid repeats or palindromes
  - Avoid high GC stretches at 5' end

### Validation Checks

The wizard automatically checks for:
- ❌ Hairpin formation (self-complementarity)
- ❌ Self-dimer formation (primer binding to itself)
- ❌ Hetero-dimer formation (primer pairs binding to each other)
- ❌ Secondary structures
- ✅ GC content within range
- ✅ GC clamp presence
- ✅ Optimal Tm range

## GUI Application

The GUI provides the most user-friendly experience:

1. **Launch**: Run `python gibson_gui.py`
2. **Input Sequences**: 
   - Paste or load vector sequence (FASTA/GenBank/text)
   - Paste or load insert sequence
3. **Set Parameters**:
   - Cut site range (e.g., 100-500 bp)
   - Homology length (15-40 nt)
   - Number of top results to show
4. **Optimize**: Click "Find Best Cut Sites"
5. **View Results**: 
   - Compare scores in results table
   - View detailed primer information
   - Export primers to file

### Features:
- 🔍 **Automatic optimization** - finds best cut sites in your specified range
- 📊 **Score-based ranking** - primers scored by Tm, GC, dimers, etc.
- 📝 **Detailed validation** - see warnings and errors for each primer
- 💾 **Export options** - save selected or all results
- 📁 **File loading** - supports FASTA, GenBank, plain text

## Command-Line Examples

### Example 1: Simple Binary Assembly

```python
# Input vector sequence (2500 bp plasmid)
Vector: ATCG... (paste your vector sequence)
Cut site: position 1250

# Input insert sequence (1000 bp gene)
Insert: GCTA... (paste your insert sequence)

# Choose homology length
Homology: 25 nt (default)
```

The wizard will generate:
- **2 Vector primers**: For linearizing the vector with insert homology
- **2 Insert primers**: For amplifying the insert with vector homology

### Example 2: Using Sequence Files

```python
from utils import read_sequence_file

# Read from FASTA or plain text file
vector_seq = read_sequence_file("vector.fasta")
insert_seq = read_sequence_file("insert.txt")
```

### Example 3: Manual Primer Design

```python
from primer_design import PrimerDesigner

designer = PrimerDesigner()

# Design a Gibson primer
template = "ATCGATCGATCG..."
homology = "GCTAGCTAGCTAGCTAGCTA"  # 20 nt homology

primer, validation = designer.design_gibson_primer(
    template=template,
    start=0,
    homology_seq=homology,
    direction='forward'
)

print(f"Primer: {primer}")
print(f"Tm: {validation['annealing_region']['tm']}°C")
print(f"Valid: {validation['valid']}")
```

## Output

### Primer Information

For each primer, the wizard provides:

```
1️⃣ Vector Forward Primer (VF):
   5'- GCTAGCTAGCTAGCTAGCTA ATCGATCGATCGATCG -3'
   Length: 36 nt
   
📊 Validation Results:
  Length: 36 nt
  GC Content: 55.6%
  Tm: 62.3°C
  GC Clamp: ✓ Yes
  
  Annealing Region:
    Length: 16 nt
    Tm: 58.4°C
    Sequence: ATCGATCGATCGATCG
  
  Homology Region:
    Length: 20 nt
    Tm: 56.2°C
    Sequence: GCTAGCTAGCTAGCTAGCTA
  
✅ All checks passed!
```

### Assembly Protocol

The wizard also provides:
- PCR conditions (annealing temperature, extension time)
- Gibson assembly protocol (Takara InFusion)
- Tips for optimization
- Option to save primers to file

## File Structure

```
gibson_assembly2.0/
│
├── primer_design.py      # Core primer design module
├── gibson_wizard.py      # Interactive wizard (main script)
├── utils.py              # Utility functions
├── README.md             # This file
└── examples/             # Example sequences (optional)
```

## Module Documentation

### `primer_design.py`

Main class: `PrimerDesigner`

Key methods:
- `calculate_tm_nearest_neighbor()`: Accurate Tm calculation
- `design_annealing_region()`: Design PCR annealing part
- `design_gibson_primer()`: Complete Gibson primer with homology
- `validate_primer()`: Comprehensive validation
- `check_hairpin()`: Hairpin structure detection
- `check_self_dimer()`: Self-dimer detection
- `check_hetero_dimer()`: Hetero-dimer detection

### `gibson_wizard.py`

Main class: `BinaryAssemblyWizard`

Interactive wizard that:
- Guides through input collection
- Designs all required primers
- Validates designs
- Generates reports
- Saves results to file

### `utils.py`

Utility functions:
- `parse_fasta()`: Parse FASTA format
- `read_sequence_file()`: Read sequences from files
- `format_sequence_block()`: Format sequences for display
- `suggest_pcr_conditions()`: Generate PCR protocols
- `generate_primer_order_sheet()`: Create order sheets

## Tips for Best Results

### Template Preparation
- Use high-quality, purified plasmid DNA
- Verify sequence by sequencing before designing primers
- Check for repeats or secondary structures in target regions

### Primer Design
- Keep total primer length ≤60 nt when possible
- For difficult templates, try different homology lengths (20-30 nt range)
- Ensure unique homology sequences to avoid mis-assembly

### PCR Optimization
- Use high-fidelity polymerase (Q5, Phusion, PrimeSTAR)
- Optimize annealing temperature (try Tm ± 3°C)
- Use 25-35 cycles to minimize errors
- Gel-purify products to remove template DNA
- Consider DpnI digest to remove methylated template

### Gibson Assembly
- Use equimolar ratios (1:1 to 1:2 vector:insert)
- Total DNA: 50-200 ng
- Takara InFusion: 50°C for 15 min
- Transform immediately after assembly

## Troubleshooting

### Low Assembly Efficiency
- Increase homology length (try 25-30 nt)
- Check primer quality and concentration
- Ensure complete removal of template DNA
- Verify PCR products are correct size and clean

### Non-specific Amplification
- Lower annealing temperature
- Reduce primer concentration
- Use touchdown PCR
- Redesign primers with higher Tm

### Secondary Structures
- Redesign primers to avoid high GC content at 5' end
- Keep homology regions simple (avoid repeats)
- Use shorter homology if needed (20 nt minimum)

## References

- Takara Bio InFusion® HD Cloning Kit User Manual
- Gibson et al. (2009) Nature Methods 6, 343-345
- SantaLucia & Hicks (2004) Annu. Rev. Biophys. Biomol. Struct. 33:415-40

## License

This tool is provided as-is for academic and research use.

## Contributing

Feel free to submit issues, suggestions, or improvements!

## Author

Created for molecular biology researchers working with Gibson assembly.

---

**Happy Cloning! 🧬🔬**
