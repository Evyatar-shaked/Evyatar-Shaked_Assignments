"""
Configuration file for Gibson Assembly Primer Designer
Modify these parameters to customize primer design behavior
"""

# PCR Primer Design Parameters
ANNEALING_LENGTH_MIN = 18  # Minimum annealing region length (nt)
ANNEALING_LENGTH_MAX = 25  # Maximum annealing region length (nt)
TARGET_TM_MIN = 55  # Minimum target Tm for annealing region (°C)
TARGET_TM_MAX = 65  # Maximum target Tm for annealing region (°C)
OPTIMAL_TM = 60  # Optimal Tm for primer design (°C)

# Gibson Homology Parameters
HOMOLOGY_LENGTH_MIN = 15  # Minimum homology region length (nt)
HOMOLOGY_LENGTH_MAX = 40  # Maximum homology region length (nt)
HOMOLOGY_LENGTH_RECOMMENDED_MIN = 20  # Recommended minimum (nt)
HOMOLOGY_LENGTH_RECOMMENDED_MAX = 30  # Recommended maximum (nt)
DEFAULT_HOMOLOGY_LENGTH = 25  # Default homology length (nt)

# Sequence Quality Parameters
GC_CONTENT_MIN = 40  # Minimum GC content (%)
GC_CONTENT_MAX = 60  # Maximum GC content (%)
MAX_TM_DIFFERENCE = 3  # Maximum Tm difference between primer pairs (°C)

# Secondary Structure Detection
MIN_HAIRPIN_STEM = 4  # Minimum stem length for hairpin detection (bp)
MIN_DIMER_MATCH = 4  # Minimum consecutive matches for dimer detection (bp)

# PCR Conditions (for protocol generation)
PRIMER_CONCENTRATION = 0.5  # Primer concentration (µM)
NA_CONCENTRATION = 50  # Na+ concentration (mM)
MG_CONCENTRATION = 1.5  # Mg2+ concentration (mM)

# Default PCR cycling conditions
PCR_DENATURATION_TEMP = 98  # Denaturation temperature (°C)
PCR_DENATURATION_TIME = 10  # Denaturation time (seconds)
PCR_ANNEALING_TIME = 30  # Annealing time (seconds)
PCR_EXTENSION_TEMP = 72  # Extension temperature (°C)
PCR_EXTENSION_RATE = 30  # Extension rate (seconds per kb)
PCR_FINAL_EXTENSION_TIME = 300  # Final extension time (seconds)
PCR_CYCLES = 30  # Number of PCR cycles

# Gibson Assembly Conditions (Takara InFusion)
ASSEMBLY_TEMPERATURE = 50  # Assembly incubation temperature (°C)
ASSEMBLY_TIME = 15  # Assembly incubation time (minutes)
ASSEMBLY_DNA_MIN = 50  # Minimum total DNA (ng)
ASSEMBLY_DNA_MAX = 200  # Maximum total DNA (ng)
ASSEMBLY_VECTOR_INSERT_RATIO = "1:1 to 1:2"  # Recommended ratio

# Primer Synthesis
DEFAULT_SYNTHESIS_SCALE = "25nm"  # Default synthesis scale
PURIFICATION_THRESHOLD = 60  # Length above which HPLC is recommended (nt)

# Validation Strictness
# Set to True for strict validation (warnings become errors)
STRICT_GC_CONTENT = False
STRICT_TM_RANGE = False
STRICT_HAIRPIN = False
STRICT_DIMER = False

# Output Formatting
SEQUENCE_DISPLAY_WIDTH = 60  # Characters per line for sequence display
SEQUENCE_GROUP_SIZE = 10  # Characters per group in sequence display

# File Output
DEFAULT_OUTPUT_FILENAME = "gibson_primers.txt"
SAVE_VALIDATION_DETAILS = True  # Include full validation in output files

# Advanced Options
ALLOW_CIRCULAR_VECTOR = True  # Allow circular vector handling
CONTEXT_LENGTH = 50  # Number of bp to show around cut site for context
OPTIMIZE_FOR_COST = False  # If True, prefer shorter primers when possible
PREFER_AT_RICH_3PRIME = False  # If True, avoid GC clamp (for AT-rich templates)

# Warning/Error Messages
VERBOSE_VALIDATION = True  # Show detailed validation messages
SHOW_WARNINGS = True  # Display warnings during design
SHOW_TIPS = True  # Display optimization tips

# Color coding for terminal output (if supported)
USE_COLOR = False  # Set to True if your terminal supports color codes

"""
Usage in code:

from config import *

designer = PrimerDesigner()
designer.min_annealing_length = ANNEALING_LENGTH_MIN
designer.max_annealing_length = ANNEALING_LENGTH_MAX
designer.target_tm_range = (TARGET_TM_MIN, TARGET_TM_MAX)
# ... etc
"""

# Presets for different protocols
PRESETS = {
    "takara_infusion": {
        "assembly_temp": 50,
        "assembly_time": 15,
        "homology_length": 25,
        "description": "Takara InFusion HD Cloning"
    },
    "neb_hifi": {
        "assembly_temp": 50,
        "assembly_time": 15,
        "homology_length": 20,
        "description": "NEB HiFi DNA Assembly"
    },
    "neb_gibson": {
        "assembly_temp": 50,
        "assembly_time": 60,
        "homology_length": 30,
        "description": "NEB Gibson Assembly Master Mix"
    },
    "custom": {
        "assembly_temp": 50,
        "assembly_time": 15,
        "homology_length": 25,
        "description": "Custom Gibson Assembly"
    }
}

# Select active preset
ACTIVE_PRESET = "takara_infusion"
