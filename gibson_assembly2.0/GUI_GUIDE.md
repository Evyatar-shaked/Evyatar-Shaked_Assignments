# 🧬 Gibson Assembly Primer Designer - GUI Quick Start

## Launch the GUI

```bash
python gibson_gui.py
```

## Step-by-Step Guide

### 1. Input Sequences Tab

**Vector Sequence:**
- Paste your vector (plasmid) sequence directly into the text box
- OR click "Load from File" to import FASTA, GenBank, or text files
- Sequence will be automatically cleaned (spaces/newlines removed)

**Insert Sequence:**
- Paste your insert (gene/fragment) sequence
- OR click "Load from File" to import

**Minimum Requirements:**
- Vector: ≥100 bp
- Insert: ≥20 bp

### 2. Set Parameters

**Cut Site Range:**
- **Start (bp)**: Beginning of range to scan (e.g., 100)
- **End (bp)**: End of range to scan (e.g., 500)
- Program will test every position in this range

**Homology Length:**
- Length of overlap regions (15-40 nt)
- Default: 25 nt (recommended for Takara InFusion)

**Top Results:**
- Number of best results to show (1-20)
- Default: 5

**Scan Step:**
- How many bp to skip between tests (1-50)
- Default: 10 bp (faster)
- Use 1 bp for most thorough search

### 3. Run Optimization

Click **"🔍 Find Best Cut Sites"**

The program will:
- Scan all positions in your range
- Design primers for each position
- Score each design based on:
  - Tm (prefer 55-65°C)
  - GC content (prefer 40-60%)
  - Hairpins and dimers (fewer is better)
  - Primer length (prefer 40-60 nt)
  - Tm compatibility between pairs

Progress bar shows activity (may take 30 seconds to few minutes)

### 4. View Results Tab

**Results Table:**
- **Rank**: Best to worst (1 = best)
- **Cut Site**: Position in vector (bp)
- **Overall Score**: Combined score (0-100, higher is better)
- **Vector Score**: Score for vector primers
- **Insert Score**: Score for insert primers
- **Warnings**: Number of warnings
- **Errors**: Number of errors

**Tips:**
- ✅ Scores >80 = Excellent
- ✅ Scores 60-80 = Good
- ⚠️ Scores 40-60 = Acceptable
- ❌ Scores <40 = Poor

Click any row to view detailed primer information

### 5. Primer Details Tab

Shows detailed information for selected result:

**For Each Primer:**
- Full sequence (5' → 3')
- Length and GC content
- Tm (annealing region)
- All warnings and errors

**Color Coding:**
- 🟢 Green = No issues
- 🟡 Orange = Warnings (review recommended)
- 🔴 Red = Errors (must fix)

### 6. Export Results

**Export Selected:**
- Save currently selected result to file
- Includes all 4 primers with details

**Export All:**
- Save all top results to one file
- Great for comparing multiple options

## Example Workflow

### Quick Test (Fast)
```
Vector: 2500 bp plasmid
Insert: 720 bp GFP
Range: 1000-1500 bp
Homology: 25 nt
Top Results: 5
Scan Step: 20 bp
```
Time: ~30 seconds

### Thorough Search (Slower but More Accurate)
```
Vector: 2500 bp plasmid
Insert: 720 bp GFP
Range: 100-2000 bp
Homology: 25 nt
Top Results: 10
Scan Step: 5 bp
```
Time: ~2-5 minutes

### Fine-tuning Best Result
Once you find a good range:
```
Range: 1200-1300 bp (narrow range)
Scan Step: 1 bp (test every position)
Top Results: 3
```
Time: ~30 seconds

## Tips for Best Results

### 1. Choose Good Range
- Avoid regions with:
  - Restriction sites you need
  - Important promoter/terminator sequences
  - Repetitive sequences
- Start with broader range (200-500 bp window)
- Refine to narrower range for fine-tuning

### 2. Optimize Parameters
- **For speed**: Step = 10-20 bp
- **For accuracy**: Step = 1-5 bp
- **Standard homology**: 25 nt
- **Difficult sequences**: Try 20 or 30 nt homology

### 3. Interpret Scores
- **Score = 90-100**: Perfect! Order immediately
- **Score = 80-90**: Excellent, ready to use
- **Score = 70-80**: Good, minor warnings acceptable
- **Score = 60-70**: Acceptable, review warnings
- **Score < 60**: Consider different range or settings

### 4. Review Top Results
- Don't just pick #1 - review top 3-5
- Check warnings for each
- Consider:
  - Position convenience (avoid features)
  - Number of warnings
  - Specific warning types

### 5. Validate Before Ordering
- Check primer sequences don't hit restriction sites
- Verify homology regions are unique
- Confirm Tm values are reasonable for your polymerase

## Troubleshooting

### "No valid cut sites found"
- **Cause**: No good primers in range
- **Solution**: 
  - Try different range
  - Adjust homology length
  - Check sequence quality

### "Very low scores (<50)"
- **Cause**: Difficult sequence region
- **Solution**:
  - Scan larger range
  - Try different homology length
  - Consider alternative vector linearization site

### GUI freezes during scan
- **Cause**: Large range + small step = many calculations
- **Solution**:
  - Use larger step size (10-20 bp)
  - Reduce range size
  - Be patient (it's working!)

### All results have many warnings
- **Cause**: Challenging sequence composition
- **Solution**:
  - Pick result with fewest warnings
  - Review specific warnings
  - Some warnings are OK (Tm slightly outside range, etc.)

## Keyboard Shortcuts

- **Ctrl+O**: Load file (when focused on text area)
- **Ctrl+A**: Select all (in text areas)
- **Ctrl+C**: Copy
- **Ctrl+V**: Paste

## Next Steps

1. **Order primers** from your preferred vendor
2. **PCR amplify**:
   - Vector primers → Linearized vector
   - Insert primers → Amplified insert
3. **Gel purify** both products
4. **Gibson assembly** following Takara InFusion protocol
5. **Transform** and select colonies

## Need Help?

- Check full documentation: `README.md`
- Run command-line examples: `python examples.py`
- View test sequences: `test_sequences.md`

---

**Ready? Launch the GUI:** `python gibson_gui.py`

🧬 Happy Cloning! 🔬
