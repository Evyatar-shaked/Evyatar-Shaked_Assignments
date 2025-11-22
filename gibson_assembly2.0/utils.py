"""
Utility functions for Gibson assembly primer design
Uses Biopython for sequence I/O
"""

from typing import List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO


def parse_fasta(fasta_text: str) -> dict:
    """
    Parse FASTA format sequence using Biopython
    
    Args:
        fasta_text: FASTA formatted text
    
    Returns:
        Dictionary with sequence name and sequence
    """
    sequences = {}
    fasta_io = StringIO(fasta_text)
    for record in SeqIO.parse(fasta_io, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    return sequences


def read_sequence_file(filename: str) -> str:
    """
    Read sequence from file using Biopython (supports FASTA, GenBank, etc.)
    
    Args:
        filename: Path to sequence file
    
    Returns:
        DNA sequence as string
    """
    try:
        # Try FASTA first
        try:
            record = SeqIO.read(filename, "fasta")
            return str(record.seq).upper()
        except:
            pass
        
        # Try GenBank
        try:
            record = SeqIO.read(filename, "genbank")
            return str(record.seq).upper()
        except:
            pass
        
        # Try plain text
        with open(filename, 'r') as f:
            content = f.read()
        sequence = ''.join(content.split()).upper()
        
        # Validate
        valid_bases = set('ATGC')
        if all(base in valid_bases for base in sequence):
            return sequence
        else:
            raise ValueError("Invalid DNA sequence in file")
    
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filename}")
    except Exception as e:
        raise Exception(f"Error reading file: {e}")


def format_sequence_block(sequence: str, width: int = 60, group: int = 10) -> str:
    """
    Format sequence in readable blocks
    
    Args:
        sequence: DNA sequence
        width: Characters per line
        group: Characters per group (separated by space)
    
    Returns:
        Formatted sequence string
    """
    lines = []
    for i in range(0, len(sequence), width):
        line_seq = sequence[i:i+width]
        # Add spaces between groups
        grouped = ' '.join([line_seq[j:j+group] for j in range(0, len(line_seq), group)])
        lines.append(f"{i+1:6d}  {grouped}")
    return '\n'.join(lines)


def calculate_product_sizes(vector_length: int, insert_length: int, 
                           cut_site: int) -> Tuple[int, int]:
    """
    Calculate expected PCR product sizes
    
    Args:
        vector_length: Length of vector in bp
        insert_length: Length of insert in bp
        cut_site: Position of cut site in vector (0-based)
    
    Returns:
        (linearized_vector_size, amplified_insert_size)
    """
    # Vector is linearized at cut site
    linearized_vector = vector_length
    
    # Insert gets the insert sequence
    amplified_insert = insert_length
    
    return linearized_vector, amplified_insert


def estimate_primer_cost(primer_length: int, scale: str = "25nm") -> float:
    """
    Estimate primer synthesis cost
    
    Args:
        primer_length: Length of primer in nucleotides
        scale: Synthesis scale (25nm, 50nm, 100nm, 200nm, 1umol)
    
    Returns:
        Estimated cost in USD (rough estimate)
    """
    base_cost = {
        "25nm": 5.0,
        "50nm": 7.0,
        "100nm": 10.0,
        "200nm": 15.0,
        "1umol": 25.0
    }
    
    cost_per_base = {
        "25nm": 0.15,
        "50nm": 0.20,
        "100nm": 0.25,
        "200nm": 0.30,
        "1umol": 0.40
    }
    
    scale = scale.lower()
    if scale not in base_cost:
        scale = "25nm"
    
    # Base cost + cost per nucleotide
    total_cost = base_cost[scale] + (primer_length * cost_per_base[scale])
    
    # Premium for long primers (>60 nt)
    if primer_length > 60:
        total_cost *= 1.5
    
    return round(total_cost, 2)


def calculate_assembly_efficiency_score(homology_length: int, gc_content: float,
                                        tm_difference: float) -> float:
    """
    Calculate predicted assembly efficiency score (0-100)
    
    Args:
        homology_length: Length of homology region
        gc_content: GC content percentage
        tm_difference: Tm difference between primer pairs
    
    Returns:
        Efficiency score (0-100)
    """
    score = 100.0
    
    # Homology length factor
    if homology_length < 15:
        score *= 0.5
    elif homology_length < 20:
        score *= 0.7
    elif homology_length > 40:
        score *= 0.9
    
    # GC content factor
    if gc_content < 30 or gc_content > 70:
        score *= 0.6
    elif gc_content < 40 or gc_content > 60:
        score *= 0.8
    
    # Tm difference factor
    if tm_difference > 5:
        score *= 0.7
    elif tm_difference > 3:
        score *= 0.85
    
    return round(score, 1)


def generate_primer_order_sheet(primers: dict, scale: str = "25nm") -> str:
    """
    Generate a formatted primer order sheet
    
    Args:
        primers: Dictionary of primer names and sequences
        scale: Synthesis scale
    
    Returns:
        Formatted order sheet as string
    """
    lines = []
    lines.append("=" * 80)
    lines.append("PRIMER ORDER SHEET")
    lines.append("=" * 80)
    lines.append(f"Scale: {scale}")
    lines.append(f"Purification: Standard desalting (or HPLC for >60 nt)")
    lines.append("")
    
    total_cost = 0
    for i, (name, sequence) in enumerate(primers.items(), 1):
        cost = estimate_primer_cost(len(sequence), scale)
        total_cost += cost
        
        lines.append(f"{i}. {name}")
        lines.append(f"   Length: {len(sequence)} nt")
        lines.append(f"   Sequence: 5'- {sequence} -3'")
        lines.append(f"   Est. Cost: ${cost:.2f}")
        lines.append("")
    
    lines.append("-" * 80)
    lines.append(f"Total Estimated Cost: ${total_cost:.2f}")
    lines.append("=" * 80)
    
    return '\n'.join(lines)


def validate_circular_vector(sequence: str) -> bool:
    """
    Check if vector sequence appears to be circular (complete plasmid)
    
    Args:
        sequence: DNA sequence
    
    Returns:
        True if likely circular, False otherwise
    """
    # Heuristic checks
    seq = sequence.upper()
    
    # Check for common plasmid features
    has_ori = 'TTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGTTTGCCGGATC' in seq  # pUC ori
    has_ampr = 'ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA' in seq  # AmpR
    
    # Check length (plasmids typically 2-15 kb)
    reasonable_length = 1000 <= len(seq) <= 20000
    
    return reasonable_length and (has_ori or has_ampr or len(seq) > 2000)


def suggest_pcr_conditions(tm: float, product_size: int) -> dict:
    """
    Suggest PCR cycling conditions based on Tm and product size
    
    Args:
        tm: Melting temperature of primers
        product_size: Expected product size in bp
    
    Returns:
        Dictionary with PCR conditions
    """
    # Annealing temperature: Tm - 5°C
    annealing_temp = max(50, tm - 5)
    
    # Extension time: 1 min per kb for high-fidelity polymerase
    extension_time = max(15, (product_size / 1000) * 60)
    
    # Cycles
    cycles = 30 if product_size < 5000 else 35
    
    return {
        'denaturation': {'temp': 98, 'time': 10},  # seconds
        'annealing': {'temp': round(annealing_temp, 1), 'time': 30},
        'extension': {'temp': 72, 'time': int(extension_time)},
        'cycles': cycles,
        'final_extension': {'temp': 72, 'time': 300}  # 5 min
    }


def format_pcr_protocol(conditions: dict, primer_names: Tuple[str, str]) -> str:
    """
    Format PCR protocol as readable text
    
    Args:
        conditions: PCR conditions dictionary
        primer_names: Tuple of (forward, reverse) primer names
    
    Returns:
        Formatted protocol string
    """
    fwd, rev = primer_names
    
    protocol = f"""
PCR PROTOCOL
============
Primers: {fwd} + {rev}

Template: 10-100 ng
Primer concentration: 0.5 µM each
High-fidelity polymerase (e.g., Q5, Phusion, PrimeSTAR)

Cycling Conditions:
-------------------
1. Initial denaturation:  {conditions['denaturation']['temp']}°C for {conditions['denaturation']['time']} sec

2. Cycles (x{conditions['cycles']}):
   - Denaturation:  {conditions['denaturation']['temp']}°C for {conditions['denaturation']['time']} sec
   - Annealing:     {conditions['annealing']['temp']}°C for {conditions['annealing']['time']} sec
   - Extension:     {conditions['extension']['temp']}°C for {conditions['extension']['time']} sec

3. Final extension:  {conditions['final_extension']['temp']}°C for {conditions['final_extension']['time']} sec

4. Hold at 4°C
"""
    return protocol
