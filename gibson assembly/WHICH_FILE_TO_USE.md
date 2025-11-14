# 📁 FILE GUIDE - What Each File Does

**You have 11 files. Don't panic! You only need to use 1-2 of them.** 

---

## ⭐ **START HERE** (The one file you need!)

### `gibson_primer_designer.ipynb` 
**👉 THIS IS THE MAIN TOOL - USE THIS ONE!**

**What it does:**
- Design Gibson assembly primers for binary cloning
- Has all functions built-in
- Includes 7 examples you can run
- Shows quality analysis (Tm, GC%, hairpins, dimers)
- Interactive - just open and run cells

**How to use:**
1. Open it in Jupyter/VS Code
2. Run cells 1-6 (loads all functions)
3. Try Examples 1-7, or customize the last cell with your sequences
4. Get your 4 primers ready to order!

**When to use:** Always! This is your main tool.

---

## 📚 **READ THESE** (To understand how it works)

### `START_HERE.md`
**Quick overview of everything**

**What it does:**
- Explains what all the files are
- Answers your questions about Tm, hairpins, dimers
- Shows what's accurate vs. simplified
- Tells you when to validate

**How to use:** Read it first to get oriented

**When to use:** When you're confused or want to understand the tool

---

### `QUALITY_CHECKS_EXPLAINED.md`
**Explains what gets checked and how accurate it is**

**What it does:**
- Shows how Tm is calculated (Nearest-Neighbor method)
- Explains hairpin/dimer detection (simplified)
- Gives interpretation guide (what's good/warning/bad)
- Compares this tool vs. professional tools

**How to use:** Read when you want to know "is this accurate enough?"

**When to use:** Before ordering primers or when you see warnings

---

### `ALGORITHM_EXPLAINED.md`
**Deep dive into how primers are designed**

**What it does:**
- Explains the construction algorithm (not search)
- Shows how annealing length is optimized
- Details the Tm calculation method
- Lists what's checked and what's not

**How to use:** Read if you're curious about the technical details

**When to use:** When you want to understand the science behind it

---

### `VISUAL_GUIDE.md`
**Diagrams and flowcharts**

**What it does:**
- Visual explanation of primer structure
- Flowcharts showing the algorithm
- Decision trees for when to use what
- ASCII art diagrams

**How to use:** Read if you're a visual learner

**When to use:** When text explanations aren't clicking

---

## 📖 **REFERENCE** (Look things up when needed)

### `README.md`
**Technical documentation**

**What it does:**
- Installation instructions
- API reference (function parameters)
- Code examples
- Feature list

**How to use:** Look up function parameters or syntax

**When to use:** When you need to remember how to call a function

---

## 🔧 **ALTERNATIVE VERSIONS** (Usually don't need these)

### `gibson_primer_design.py`
**Python module version (same as notebook but as .py file)**

**What it does:**
- Same functionality as the notebook
- Can import into other Python scripts
- For integration into workflows

**How to use:** 
```python
from gibson_primer_design import GibsonPrimerDesigner
```

**When to use:** 
- If you want to import it into another script
- If you prefer .py files over notebooks
- **Most people won't need this** - just use the notebook!

---

### `gibson_demo.ipynb`
**Simpler demo notebook**

**What it does:**
- Shorter version with basic examples
- Imports from `gibson_primer_design.py`
- Less comprehensive than `gibson_primer_designer.ipynb`

**When to use:**
- **You probably don't need this** - the main notebook is better!
- Only useful if you want a minimal example

---

## 🧪 **TESTING/SETUP** (Background files)

### `test_gibson.py`
**Automated tests**

**What it does:**
- Tests the code to make sure it works
- Runs 4 test scenarios
- Validates functionality

**How to use:**
```bash
python test_gibson.py
```

**When to use:** 
- After installation to verify everything works
- **Most people don't need to run this**

---

### `requirements.txt`
**Package list**

**What it does:**
- Lists required Python packages (pydna, biopython)

**How to use:**
```bash
pip install -r requirements.txt
```

**When to use:** 
- First time setup
- To install dependencies

---

## 🗑️ **CAN PROBABLY IGNORE**

### `gibson_assembly_planner.ipynb` (your original file)
**Your original notebook before I created the new one**

**What it does:**
- Had basic primer design
- No quality checks

**When to use:** You can probably delete this - the new notebook is better!

---

### `gibson_assemby.ipynb` (typo in name)
**Another notebook (possibly empty or old version)**

**When to use:** Check if it has anything important, otherwise can delete

---

## 🎯 **SIMPLE RECOMMENDATION**

### **For 99% of people, just use these:**

1. **`gibson_primer_designer.ipynb`** ← Do your work here
2. **`QUALITY_CHECKS_EXPLAINED.md`** ← Read to understand results

### **That's it!** 

---

## 📋 **QUICK START (3 steps)**

```
Step 1: Install packages
  → pip install pydna biopython

Step 2: Open the main notebook
  → gibson_primer_designer.ipynb

Step 3: Run cells 1-6, then try examples
  → Done! Design your primers!
```

---

## 🗂️ **FILE ORGANIZATION**

```
YOUR FOLDER:
│
├─ 🟢 USE THIS: gibson_primer_designer.ipynb (MAIN TOOL)
│
├─ 📖 READ THESE (if you want to understand):
│   ├─ START_HERE.md
│   ├─ QUALITY_CHECKS_EXPLAINED.md
│   ├─ ALGORITHM_EXPLAINED.md
│   └─ VISUAL_GUIDE.md
│
├─ 📚 REFERENCE (lookup when needed):
│   └─ README.md
│
├─ 🔧 ALTERNATIVE (usually don't need):
│   ├─ gibson_primer_design.py
│   └─ gibson_demo.ipynb
│
├─ ⚙️ SETUP/TESTING (background):
│   ├─ requirements.txt
│   └─ test_gibson.py
│
└─ 🗑️ OLD FILES (probably can ignore):
    ├─ gibson_assembly_planner.ipynb
    └─ gibson_assemby.ipynb
```

---

## ❓ **WHICH FILE FOR WHICH TASK?**

### "I want to design primers"
→ **`gibson_primer_designer.ipynb`**

### "I want to understand how it works"
→ **`QUALITY_CHECKS_EXPLAINED.md`** or **`ALGORITHM_EXPLAINED.md`**

### "I see warnings, are they serious?"
→ **`QUALITY_CHECKS_EXPLAINED.md`** (has interpretation guide)

### "How do I call this function?"
→ **`README.md`** (has API reference)

### "I want to import into my own script"
→ **`gibson_primer_design.py`**

### "I want to install it"
→ **`requirements.txt`** (run `pip install -r requirements.txt`)

---

## 💡 **BOTTOM LINE**

**Too many files? Just focus on this:**

```
┌─────────────────────────────────────────┐
│                                         │
│  Open: gibson_primer_designer.ipynb    │
│                                         │
│  Read: QUALITY_CHECKS_EXPLAINED.md     │
│        (if you want)                    │
│                                         │
│  Ignore: Everything else for now!      │
│                                         │
└─────────────────────────────────────────┘
```

**That's all you need to get started!** 🎉

The other files are just documentation and alternatives. You can explore them later if you want, but they're not necessary to use the tool.

---

## 🤔 **CONFUSED? HERE'S THE MINIMUM:**

1. Open `gibson_primer_designer.ipynb`
2. Run the cells
3. Put in your sequences
4. Get your primers
5. Done!

**Seriously, it's that simple.** All the other files are just extras to help you understand or use it in different ways.
