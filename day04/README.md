# DNA BLAST Search Tool

A Python application that performs DNA sequence alignment using NCBI BLAST API with a user-friendly graphical interface.

## Features

- 🧬 Search DNA sequences against NCBI's nucleotide database
- 🔍 Filter results by organism (dropdown with common organisms)
- 📊 Display detailed alignment statistics (identity %, E-value, scores)
- 🖥️ Easy-to-use GUI built with Tkinter
- 📝 View alignment details including matches and gaps
- ⚡ Asynchronous search to keep UI responsive

## Requirements

- Python 3.7 or higher
- Biopython library
- Tkinter (included with Python on Windows)

## Installation

1. **Clone or download this repository**

2. **Install dependencies:**

   ```powershell
   pip install -r requirements.txt
   ```

   Or install Biopython directly:
   
   ```powershell
   pip install biopython
   ```

## Usage

### Running the GUI Application

Launch the graphical interface:

```powershell
python blast_gui.py
```

### Using the GUI:

1. **Enter DNA Sequence**: Paste or type your DNA sequence in the text area
2. **Select Organism**: Choose an organism from the dropdown menu or select "All organisms"
3. **Set Max Results**: Choose how many results to return (1-50)
4. **Click "Search BLAST"**: Wait for results (typically 1-2 minutes)
5. **View Results**: Scroll through the detailed alignment results

### Example Sequence

Click the "Load Example Sequence" button to load a fragment of the human insulin gene for testing.

### Command Line Usage

You can also use the BLAST search module programmatically:

```python
from blast_search import BLASTSearcher

# Create searcher instance
searcher = BLASTSearcher()

# Define your DNA sequence
sequence = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTC"

# Validate sequence
is_valid, validated_seq = searcher.validate_sequence(sequence)

if is_valid:
    # Perform search
    results = searcher.search_sequence(
        validated_seq, 
        organism="Homo sapiens",
        max_results=10
    )
    
    # Display results
    for result in results:
        print(searcher.format_result_summary(result))
```

## How It Works

1. **Input Validation**: Checks that the DNA sequence contains only valid nucleotides (A, T, C, G, N)
2. **BLAST Query**: Submits the sequence to NCBI's BLAST service via Biopython's NCBIWWW module
3. **XML Parsing**: Parses the XML results returned by NCBI
4. **Result Display**: Shows alignment statistics including:
   - Sequence title and accession number
   - Identity percentage
   - E-value (statistical significance)
   - Alignment score
   - Visual alignment preview

## Supported Organisms

The dropdown includes common model organisms:
- Homo sapiens (Human)
- Mus musculus (Mouse)
- Rattus norvegicus (Rat)
- Drosophila melanogaster (Fruit fly)
- Caenorhabditis elegans (Roundworm)
- Saccharomyces cerevisiae (Yeast)
- Escherichia coli (Bacteria)
- Arabidopsis thaliana (Plant)
- Danio rerio (Zebrafish)
- Xenopus laevis (African clawed frog)

You can also type any organism name directly into the dropdown.

## Notes

- **Search Time**: BLAST searches typically take 1-2 minutes depending on sequence length and NCBI server load
- **Internet Required**: This tool requires an active internet connection to access NCBI services
- **Rate Limiting**: NCBI has usage limits; avoid submitting too many queries in rapid succession
- **Sequence Length**: Very short sequences (< 10 bp) are not accepted

## Troubleshooting

**Error: "BLAST search failed"**
- Check your internet connection
- Verify the sequence contains only valid DNA characters
- Try again in a few minutes (NCBI servers may be busy)

**No results found**
- Try searching without organism filter
- Verify your sequence is correct
- The sequence may be too short or not have close matches in the database

**GUI not opening**
- Ensure Tkinter is installed (should be included with Python on Windows)
- Try running: `python -m tkinter` to test Tkinter installation

## File Structure

```
day04/
├── blast_search.py      # Core BLAST search functionality
├── blast_gui.py         # GUI application
├── requirements.txt     # Python dependencies
└── README.md           # This file
```

## License

This is an educational project for learning bioinformatics programming.

## Credits

- Uses [Biopython](https://biopython.org/) for NCBI BLAST integration
- Queries [NCBI BLAST](https://blast.ncbi.nlm.nih.gov/) web services
