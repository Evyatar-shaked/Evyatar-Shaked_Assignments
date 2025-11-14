# Gibson Assembly Primer Design Module

A Python module for designing primers for **binary Gibson assembly** using pydna.

## Features

- 🧬 Automatic primer design for binary Gibson assembly (one insert per step)
- 🌡️ Optimal melting temperature (Tm) calculation
- 🔄 PCR simulation and assembly verification using pydna
- 📊 Multi-step sequential assembly planning
- 🎯 Customizable homology regions and insertion sites
- ⚗️ Works with circular plasmids

## Installation

```bash
pip install pydna biopython
```

## Quick Start

```python
from gibson_primer_design import GibsonPrimerDesigner
from pydna.dseqrecord import Dseqrecord

# Create your sequences
vector = Dseqrecord("ATGC..." * 100, circular=True, name="pVector")
insert = Dseqrecord("GGGATTTCCC...", name="GeneInsert")

# Design primers
designer = GibsonPrimerDesigner(homology_length=40, target_tm=60.0)
result = designer.design_primers(vector, insert, insert_site=100)

# Print summary
designer.print_primer_summary(result)

# Get primer sequences for ordering
primers = result['primers']
print(primers['insert_fwd'])
print(primers['insert_rev'])
print(primers['vector_fwd'])
print(primers['vector_rev'])
```

## Binary Gibson Assembly Concept

In binary Gibson assembly, each step adds **one insert to a vector**. This requires **4 primers**:

1. **Insert Forward**: Amplifies insert with vector left homology tail
2. **Insert Reverse**: Amplifies insert with vector right homology tail
3. **Vector Forward**: Amplifies vector with insert left homology tail
4. **Vector Reverse**: Amplifies vector with insert right homology tail

After PCR:
- Insert PCR product has vector-compatible ends
- Vector PCR product has insert-compatible ends
- These assemble together via Gibson assembly

## API Reference

### `GibsonPrimerDesigner`

Main class for primer design.

**Parameters:**
- `homology_length` (int): Length of overlap for Gibson assembly (default: 40 bp)
- `min_anneal_length` (int): Minimum annealing region (default: 20 bp)
- `target_tm` (float): Target melting temperature (default: 60.0°C)
- `tm_tolerance` (float): Acceptable Tm deviation (default: ±5.0°C)

**Methods:**

#### `design_primers(vector, insert, insert_site)`

Design 4 primers for a single binary Gibson assembly step.

**Parameters:**
- `vector` (Dseqrecord): Vector backbone (can be circular)
- `insert` (Dseqrecord): Insert sequence to add
- `insert_site` (int): Position in vector for insertion (0-based)

**Returns:**
Dictionary containing:
- `primers`: Dict with primer sequences
- `primer_objects`: pydna Primer objects
- `pcr_products`: Simulated PCR products
- `assembly_result`: Predicted assembly product
- `primer_details`: Detailed info (Tm, lengths, etc.)

#### `plan_multi_step_assembly(initial_vector, inserts, insert_site)`

Plan sequential multi-step assembly.

**Parameters:**
- `initial_vector` (Dseqrecord): Starting vector
- `inserts` (List[Dseqrecord]): List of inserts to add sequentially
- `insert_site` (int): Position for insertions

**Returns:**
List of dictionaries, one per step, each containing primer designs.

#### `print_primer_summary(result)`

Print formatted summary of primer design results.

## Multi-Step Assembly Example

```python
# Create multiple inserts
inserts = [
    Dseqrecord("GGGAAATTT...", name="Insert1"),
    Dseqrecord("TTTGGGCCC...", name="Insert2"),
    Dseqrecord("AAACCCTTT...", name="Insert3"),
]

# Plan all steps
steps = designer.plan_multi_step_assembly(
    initial_vector=vector,
    inserts=inserts,
    insert_site=100
)

# Access primers for each step
for step in steps:
    print(f"Step {step['step']} primers:")
    for name, seq in step['primers'].items():
        print(f"  {name}: {seq}")
```

## Working with Files

Load sequences from GenBank or FASTA files:

```python
from pydna.readers import read

# Load from file
vector = read("plasmid.gb")[0]
insert = read("gene.fasta")[0]

# Design primers
result = designer.design_primers(vector, insert, insert_site=500)
```

## Customization

### Adjust Homology Length

```python
# Shorter homology (faster but less efficient)
designer = GibsonPrimerDesigner(homology_length=20)

# Longer homology (more efficient but larger primers)
designer = GibsonPrimerDesigner(homology_length=60)
```

### Control Primer Tm

```python
# Higher Tm for difficult templates
designer = GibsonPrimerDesigner(target_tm=65.0, tm_tolerance=3.0)
```

### Specify Insertion Region

```python
# The insert_site parameter controls where homology regions are extracted
result = designer.design_primers(
    vector=vector,
    insert=insert,
    insert_site=250  # Homology will be extracted around position 250
)
```

## Output Details

Each primer has detailed information:

```python
primer_info = result['primer_details']['insert_fwd']
print(f"Sequence: {primer_info['sequence']}")
print(f"Total length: {primer_info['length']} bp")
print(f"Homology tail: {primer_info['homology_length']} bp")
print(f"Annealing region: {primer_info['anneal_length']} bp")
print(f"Annealing Tm: {primer_info['anneal_tm']}°C")
print(f"Full primer Tm: {primer_info['full_tm']}°C")
```

## Tips

1. **Homology Length**: 40 bp is optimal. Range: 20-80 bp
2. **Circular Plasmids**: Module handles circular DNA correctly
3. **Tm Matching**: All primers designed for similar Tm (~60°C)
4. **Assembly Verification**: Automatically simulates PCR and assembly
5. **Multi-Step**: Each step updates vector for next insertion

## Example Output

```
================================================================================
GIBSON ASSEMBLY PRIMER DESIGN SUMMARY
================================================================================

INSERT FWD:
  Sequence (5'→3'): ATGCGTACGTTAGCCTAGGCTAGCTAGGCTTACGATGCGGGAAAATTTCCCGGGTTTAAA
  Total Length: 60 bp
  Homology Tail: 40 bp (for Gibson assembly)
  Annealing Region: 20 bp (Tm: 58.3°C)
  Full Primer Tm: 72.1°C

INSERT REV:
  Sequence (5'→3'): TAGCTAGCCTAGCTAGCCTAGCTACGTACGCATTTTAAACCGGGAAA
  Total Length: 48 bp
  Homology Tail: 40 bp (for Gibson assembly)
  Annealing Region: 8 bp (Tm: 59.8°C)
  Full Primer Tm: 68.9°C

[... vector primers ...]

================================================================================
PCR PRODUCTS:
================================================================================
Insert PCR product: 88 bp
Vector PCR product: 720 bp

✓ Assembly Product: 808 bp
```

## See Also

- [pydna documentation](https://github.com/BjornFJohansson/pydna)
- [Gibson Assembly Protocol](https://www.neb.com/protocols/2012/12/11/gibson-assembly-protocol)

## License

MIT
