# 🚀 START HERE - Gibson Assembly Primer Designer

## What is This Tool?

A **smart Gibson assembly primer designer** that automatically finds the **best cut site** in your vector for optimal primer design.

## 🎯 Two Ways to Use

### 🖥️ Option 1: GUI Application (RECOMMENDED)

```bash
pip install biopython
python gibson_gui.py
```

**Perfect for:**
- Visual learners
- First-time users
- Comparing multiple options
- File loading
- Optimization

**What it does:**
1. Load your vector and insert sequences
2. Set a range to scan (e.g., 100-500 bp)
3. Click "Find Best Cut Sites"
4. See top results ranked by score
5. Export primers

### 💻 Option 2: Command-Line Wizard

```bash
python gibson_wizard.py
```

**Perfect for:**
- Terminal users
- Scripting
- Quick designs
- Known cut sites

## 📖 Which Guide to Read?

**Just want to start?**
→ Read `GUI_GUIDE.md` (5 minutes)

**Need quick reference?**
→ Read `QUICKSTART.md` (command-line)

**Want all details?**
→ Read `README.md` (complete documentation)

**Need to find a file?**
→ Read `INDEX.md` (file navigator)

## 🎓 Quick Tutorial

### Scenario: Clone GFP into Expression Vector

1. **Launch GUI**: `python gibson_gui.py`

2. **Input Tab**:
   - Vector: Paste your 3000 bp plasmid
   - Insert: Paste your 720 bp GFP gene

3. **Set Parameters**:
   - Cut site range: 1000-1500 bp
   - Homology: 25 nt
   - Top results: 5
   - Scan step: 10 bp

4. **Click**: "🔍 Find Best Cut Sites"

5. **Results Tab**:
   - See 5 best options ranked by score
   - Click to view details

6. **Details Tab**:
   - Review all 4 primers
   - Check warnings
   - Verify Tm values

7. **Export**: Save primers to file

**Time: ~30 seconds**

## 🏆 Key Features

### What Makes This Tool Special?

✅ **Automatic Optimization**
- Scans your chosen range
- Tests every position
- Scores based on multiple criteria
- Finds THE best spot

✅ **Smart Scoring**
- Tm optimization (55-65°C)
- GC content (40-60%)
- No hairpins or dimers
- Pair compatibility
- **Score: 0-100 (higher = better)**

✅ **Multiple Results**
- See top 5 (or more) options
- Compare scores
- Review warnings
- Choose best for your needs

✅ **Biopython Powered**
- Accurate Tm calculations
- Industry-standard algorithms
- Multiple file format support

✅ **Easy Export**
- Save primers to file
- Order directly from output
- Includes validation details

## 📊 Score Interpretation

| Score | Meaning | Action |
|-------|---------|--------|
| 90-100 | Perfect! | Order immediately ✅ |
| 80-90 | Excellent | Ready to use ✅ |
| 70-80 | Good | Minor warnings OK ⚠️ |
| 60-70 | Acceptable | Review warnings ⚠️ |
| <60 | Poor | Try different range ❌ |

## 🔧 Installation

```bash
# Install Biopython
pip install biopython

# That's it! Ready to use.
```

## 📁 Project Structure

```
gibson_assembly2.0/
│
├── gibson_gui.py         ← START HERE (GUI)
├── GUI_GUIDE.md          ← How to use GUI
│
├── gibson_wizard.py      ← Command-line version
├── QUICKSTART.md         ← CLI quick start
│
├── README.md             ← Full documentation
├── INDEX.md              ← File navigator
│
├── examples.py           ← Run to see examples
├── test_design.py        ← Test the tool
│
├── primer_design.py      ← Core engine
├── optimizer.py          ← Optimization algorithm
├── utils.py              ← Helper functions
├── config.py             ← Settings
│
└── requirements.txt      ← Dependencies
```

## 🎬 Quick Start (30 seconds)

```bash
# 1. Install
pip install biopython

# 2. Launch GUI
python gibson_gui.py

# 3. Load sequences
#    - Click "Load from File" or paste directly
#    - Enter cut site range (e.g., 100-500)
#    - Click "Find Best Cut Sites"

# 4. Done!
#    - View results sorted by score
#    - Export best primers
```

## 💡 Pro Tips

1. **Use Optimization**: Don't guess cut sites - let the tool find the best one!

2. **Start Broad**: Begin with 200-500 bp range, then refine

3. **Review Top 3-5**: Compare multiple options before choosing

4. **Check Warnings**: Some are OK (minor Tm differences), others need attention

5. **Export Results**: Save for future reference

## ❓ Common Questions

**Q: How long does it take?**
A: 30 seconds to 2 minutes depending on range and step size

**Q: What if scores are low?**
A: Try different range, adjust homology length, or accept best available

**Q: Can I use known cut site?**
A: Yes! Use command-line wizard (`gibson_wizard.py`) for manual cut site

**Q: What files can I load?**
A: FASTA, GenBank, or plain text

**Q: Do I need to install anything?**
A: Just Biopython: `pip install biopython`

## 🆘 Need Help?

1. **GUI Issues**: Read `GUI_GUIDE.md`
2. **CLI Usage**: Read `QUICKSTART.md`
3. **Technical Details**: Read `README.md`
4. **Can't Find File**: Read `INDEX.md`
5. **See Examples**: Run `python examples.py`

## 🎉 Ready to Start?

### For GUI (Recommended):
```bash
python gibson_gui.py
```

### For Command-Line:
```bash
python gibson_wizard.py
```

### For Examples:
```bash
python examples.py
```

---

## 📝 Example Output

```
Results Table:
┌──────┬──────────┬─────────────┬──────────────┬───────────────┬──────────┬────────┐
│ Rank │ Cut Site │ Score       │ Vector Score │ Insert Score  │ Warnings │ Errors │
├──────┼──────────┼─────────────┼──────────────┼───────────────┼──────────┼────────┤
│  1   │   1245   │    92.5     │     94.1     │     90.9      │    2     │   0    │
│  2   │   1255   │    91.8     │     93.2     │     90.4      │    3     │   0    │
│  3   │   1235   │    90.2     │     91.8     │     88.6      │    4     │   0    │
│  4   │   1265   │    88.7     │     90.5     │     86.9      │    4     │   0    │
│  5   │   1225   │    87.3     │     89.2     │     85.4      │    5     │   0    │
└──────┴──────────┴─────────────┴──────────────┴───────────────┴──────────┴────────┘

Best Result (Rank #1):
Vector Forward:  5'- ATGGTGAGCAAGGGCGAGGAGCTG...TGCTAGCGCTATATGCGTTG -3'
Vector Reverse:  5'- CTCGGCGCGGGTCTTGTAGTTGCC...GCCAAAGCGGTCGGACAG -3'
Insert Forward:  5'- GCTGCTAGCGCTATATGCGTTGAT...ATGGTGAGCAAGGGCGAGG -3'
Insert Reverse:  5'- CGAGAACGGGTGCGCATAGAAATT...GCTCGGCGCGGGTCTTGTA -3'

Score: 92.5/100 ✅ Excellent!
```

---

**🧬 Happy Cloning! 🔬**

**Questions?** Check the documentation files or run `python examples.py` for working examples.
