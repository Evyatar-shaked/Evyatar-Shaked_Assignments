# Quick Start Guide - Gibson Assembly Primer Wizard

## Installation

1. Install Biopython:
```bash
pip install biopython
```

2. Run the Wizard

```bash
python gibson_wizard.py
```

## Step-by-Step Walkthrough

### Step 1: Enter Vector Sequence
Paste your plasmid/vector sequence when prompted:
```
Paste vector sequence: ATCGATCG...
✓ Vector sequence loaded: 2686 bp
```

### Step 2: Specify Cut Site
Enter where to linearize the vector (1-based position):
```
Cut site position (1-2686): 1343
✓ Cut site set at position 1343
```

### Step 3: Enter Insert Sequence
Paste your gene/fragment to be cloned:
```
Paste insert sequence: ATGGTGAGC...
✓ Insert sequence loaded: 720 bp
```

### Step 4: Set Homology Length
Choose homology region length (or use default):
```
Customize homology length? (y/n): n
Using default: 25 nt
```

### Step 5: Review Primers
The wizard will design 4 primers:
- **Vector Forward (VF)**: Linearizes vector, adds insert start homology
- **Vector Reverse (VR)**: Linearizes vector, adds insert end homology
- **Insert Forward (IF)**: Amplifies insert, adds vector upstream homology
- **Insert Reverse (IR)**: Amplifies insert, adds vector downstream homology

### Step 6: Save Results
```
Save primers to file? (y/n): y
Enter filename: my_primers.txt
✓ Primers saved to my_primers.txt
```

## What You Get

For each primer:
- ✅ Full sequence (5' → 3')
- ✅ Length and Tm
- ✅ GC content
- ✅ Separate annealing and homology regions
- ✅ Validation warnings (if any)

## PCR & Assembly Protocol

### PCR (2 separate reactions)

**Reaction 1: Vector Linearization**
```
Template: Vector plasmid (10-100 ng)
Primers: VF + VR (0.5 µM each)
Polymerase: High-fidelity (Q5, Phusion, PrimeSTAR)
Expected product: ~Vector length
```

**Reaction 2: Insert Amplification**
```
Template: Insert template (10-100 ng)
Primers: IF + IR (0.5 µM each)
Polymerase: High-fidelity
Expected product: ~Insert length
```

### Gibson Assembly (Takara InFusion)

```
1. Gel-purify both PCR products
2. Mix equimolar amounts (1:1 to 1:2 ratio vector:insert)
3. Total DNA: 50-200 ng
4. Add InFusion enzyme mix
5. Incubate at 50°C for 15 minutes
6. Transform 2-5 µL into competent cells
```

## Tips for Success

### Before PCR
- Verify sequences are correct
- Check primers for secondary structures
- Order primers with standard desalting (HPLC for >60 nt)

### During PCR
- Use high-fidelity polymerase (essential!)
- Optimize annealing temp: Start at Tm - 5°C
- Run gel to confirm single band of correct size
- Gel-purify to remove template DNA

### Gibson Assembly
- Use fresh, purified PCR products
- Equimolar ratios work best
- Don't exceed 200 ng total DNA
- Transform immediately after assembly
- Plate on selective media

### Troubleshooting

**No colonies?**
- Check DNA concentration
- Verify competent cell quality
- Try increasing insert:vector ratio to 2:1
- Ensure complete template removal (DpnI digest)

**Wrong insert?**
- Verify by colony PCR before sequencing
- Check homology regions are unique
- Ensure no cross-homology between fragments

**Multiple bands in PCR?**
- Lower annealing temperature
- Redesign primers with higher Tm
- Use touchdown PCR
- Check template quality

## Examples to Try

Run the examples script:
```bash
python examples.py
```

This demonstrates:
- Basic primer validation
- Gibson primer design
- Primer pair compatibility checks
- Batch validation
- Tm calculation methods

## Need Help?

Check the full README.md for:
- Detailed design rules
- API documentation
- Advanced usage
- Troubleshooting guide

## Quick Reference

| Parameter | Recommended | Acceptable Range |
|-----------|-------------|------------------|
| Annealing length | 20 nt | 18-25 nt |
| Homology length | 25 nt | 15-40 nt |
| Annealing Tm | 60°C | 55-65°C |
| GC content | 50% | 40-60% |
| Tm difference | <1°C | <3°C |

---

**Ready to start? Run:** `python gibson_wizard.py`
