"""
Gibson Assembly Primer Design Module
Handles primer design for binary Gibson assembly following Takara InFusion protocol
Uses Biopython for sequence operations and Tm calculations
"""

from typing import Dict, List, Tuple, Optional
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt


def generate_final_construct(vector_seq: str, insert_seq: str, vector_upstream_pos: int,
                            insert_fwd_pos: int = 0, insert_rev_pos: int = None) -> str:
    """
    Generate the final assembled construct sequence
    
    Args:
        vector_seq: Vector sequence
        insert_seq: Insert sequence
        vector_upstream_pos: Position of vector forward primer (linearization point)
        insert_fwd_pos: Start position of insert region to amplify
        insert_rev_pos: End position of insert region to amplify
    
    Returns:
        Final circular construct sequence
    """
    if insert_rev_pos is None:
        insert_rev_pos = len(insert_seq)
    
    # Extract insert region
    insert_region = insert_seq[insert_fwd_pos:insert_rev_pos]
    
    # Split vector at linearization point
    vector_upstream = vector_seq[:vector_upstream_pos]
    vector_downstream = vector_seq[vector_upstream_pos:]
    
    # Assemble: upstream vector + insert + downstream vector
    final_construct = vector_upstream + insert_region + vector_downstream
    
    return final_construct


class PrimerDesigner:
    """Design Gibson assembly primers for binary assembly"""
    
    def __init__(self):
        self.min_annealing_length = 18
        self.max_annealing_length = 25
        self.target_tm_range = (55, 65)
        self.min_homology_length = 15
        self.max_homology_length = 40
        self.recommended_homology_length = (20, 30)
        self.gc_content_range = (40, 60)
        self.max_tm_difference = 3
    
    def calculate_tm_nearest_neighbor(self, seq: str, primer_conc: float = 0.5, 
                                     na_conc: float = 50, mg_conc: float = 1.5) -> float:
        """
        Calculate Tm using Biopython's nearest-neighbor method
        
        Args:
            seq: DNA sequence
            primer_conc: Primer concentration in µM (default 0.5)
            na_conc: Na+ concentration in mM (default 50)
            mg_conc: Mg2+ concentration in mM (default 1.5)
        
        Returns:
            Melting temperature in °C
        """
        if len(seq) < 2:
            return 0
        
        bio_seq = Seq(seq)
        # Use Biopython's Tm_NN method with salt correction
        # Biopython expects: dnac1 in nM, Na in mM, Mg in mM
        tm = mt.Tm_NN(
            bio_seq,
            dnac1=primer_conc * 1000,  # Convert µM to nM
            Na=na_conc,  # Na+ in mM
            Mg=mg_conc,  # Mg2+ in mM
            saltcorr=7  # Owczarzy et al. 2008 correction
        )
        return round(tm, 1)
    
    def calculate_gc_content(self, seq: str) -> float:
        """Calculate GC content as percentage using Biopython"""
        if len(seq) == 0:
            return 0
        bio_seq = Seq(seq)
        return round(gc_fraction(bio_seq) * 100, 1)
    
    def has_gc_clamp(self, seq: str) -> bool:
        """Check if sequence has GC clamp (1-2 G/C at 3' end)"""
        seq = seq.upper()
        if len(seq) < 2:
            return False
        # Check last 2 bases
        last_two = seq[-2:]
        gc_count = last_two.count('G') + last_two.count('C')
        return 1 <= gc_count <= 2
    
    def check_hairpin(self, seq: str, min_stem: int = 4) -> Tuple[bool, float]:
        """
        Check for hairpin formation (sequence folding back on itself)
        Returns (has_hairpin, max_complementarity_score)
        """
        seq = seq.upper()
        rev_comp = self.reverse_complement(seq)
        
        max_score = 0
        has_significant_hairpin = False
        
        # Check for complementarity between sequence and its reverse complement
        for i in range(len(seq) - min_stem + 1):
            for j in range(len(rev_comp) - min_stem + 1):
                score = 0
                for k in range(min(len(seq) - i, len(rev_comp) - j)):
                    if seq[i + k] == rev_comp[j + k]:
                        score += 1
                    else:
                        if score >= min_stem:
                            break
                        score = 0
                    
                    if score >= min_stem:
                        has_significant_hairpin = True
                        max_score = max(max_score, score)
                        break
        
        return has_significant_hairpin, max_score
    
    def check_self_dimer(self, seq: str, min_match: int = 4) -> Tuple[bool, int]:
        """
        Check if primer can bind to itself (self-dimer)
        Returns (has_self_dimer, max_consecutive_matches)
        """
        seq = seq.upper()
        rev_seq = seq[::-1]
        
        max_matches = 0
        has_dimer = False
        
        # Slide the sequence against itself (reversed)
        for offset in range(1, len(seq)):
            consecutive = 0
            for i in range(len(seq) - offset):
                if seq[i] == rev_seq[i + offset]:
                    consecutive += 1
                    if consecutive >= min_match:
                        has_dimer = True
                        max_matches = max(max_matches, consecutive)
                else:
                    consecutive = 0
        
        return has_dimer, max_matches
    
    def check_hetero_dimer(self, seq1: str, seq2: str, min_match: int = 4) -> Tuple[bool, int]:
        """
        Check if two primers can bind to each other (hetero-dimer)
        Returns (has_hetero_dimer, max_consecutive_matches)
        """
        seq1 = seq1.upper()
        seq2_rev_comp = self.reverse_complement(seq2)
        
        max_matches = 0
        has_dimer = False
        
        # Slide seq1 against reverse complement of seq2
        for offset in range(-(len(seq2_rev_comp)-1), len(seq1)):
            consecutive = 0
            for i in range(max(0, -offset), min(len(seq1), len(seq2_rev_comp) - offset)):
                if seq1[i] == seq2_rev_comp[i + offset]:
                    consecutive += 1
                    if consecutive >= min_match:
                        has_dimer = True
                        max_matches = max(max_matches, consecutive)
                else:
                    consecutive = 0
        
        return has_dimer, max_matches
    
    def reverse_complement(self, seq: str) -> str:
        """Return reverse complement of DNA sequence using Biopython"""
        bio_seq = Seq(seq)
        return str(bio_seq.reverse_complement())
    
    def validate_primer(self, primer: str, partner_primer: Optional[str] = None) -> Dict:
        """
        Validate a primer against all design rules
        
        Returns a dictionary with validation results
        """
        results = {
            'valid': True,
            'warnings': [],
            'errors': [],
            'metrics': {}
        }
        
        primer = primer.upper()
        
        # Calculate metrics
        length = len(primer)
        gc_content = self.calculate_gc_content(primer)
        tm = self.calculate_tm_nearest_neighbor(primer)
        has_gc_clamp = self.has_gc_clamp(primer)
        has_hairpin, hairpin_score = self.check_hairpin(primer)
        has_self_dimer, self_dimer_score = self.check_self_dimer(primer)
        
        results['metrics'] = {
            'length': length,
            'gc_content': gc_content,
            'tm': tm,
            'has_gc_clamp': has_gc_clamp,
            'hairpin_score': hairpin_score,
            'self_dimer_score': self_dimer_score
        }
        
        # Check GC content
        if gc_content < self.gc_content_range[0] or gc_content > self.gc_content_range[1]:
            results['warnings'].append(
                f"GC content {gc_content}% is outside recommended range {self.gc_content_range[0]}-{self.gc_content_range[1]}%"
            )
        
        # Check GC clamp
        if not has_gc_clamp:
            results['warnings'].append("No GC clamp at 3' end (recommend 1-2 G/C bases)")
        
        # Check hairpin
        if has_hairpin:
            results['warnings'].append(f"Potential hairpin formation detected (complementarity: {hairpin_score} bp)")
        
        # Check self-dimer
        if has_self_dimer:
            results['warnings'].append(f"Potential self-dimer detected (consecutive matches: {self_dimer_score} bp)")
        
        # Check hetero-dimer if partner provided
        if partner_primer:
            has_hetero, hetero_score = self.check_hetero_dimer(primer, partner_primer)
            results['metrics']['hetero_dimer_score'] = hetero_score
            if has_hetero:
                results['warnings'].append(
                    f"Potential hetero-dimer with partner primer (consecutive matches: {hetero_score} bp)"
                )
        
        # Check Tm range
        if tm < self.target_tm_range[0] or tm > self.target_tm_range[1]:
            results['warnings'].append(
                f"Tm {tm}°C is outside recommended range {self.target_tm_range[0]}-{self.target_tm_range[1]}°C"
            )
        
        # Set overall validity
        if results['errors']:
            results['valid'] = False
        
        return results
    
    def design_annealing_region(self, template: str, start: int, direction: str = 'forward',
                               target_length: int = 20) -> Tuple[str, Dict]:
        """
        Design the annealing region of a primer
        
        Args:
            template: DNA template sequence
            start: Starting position on template
            direction: 'forward' or 'reverse'
            target_length: Target length for annealing region
        
        Returns:
            (annealing_sequence, validation_results)
        """
        template = template.upper()
        
        if direction == 'forward':
            # Extend forward from start position
            end = min(start + target_length, len(template))
            annealing_seq = template[start:end]
        else:
            # Extend backward from start position, then reverse complement
            begin = max(0, start - target_length + 1)
            annealing_seq = self.reverse_complement(template[begin:start + 1])
        
        # Try to optimize length for Tm
        best_seq = annealing_seq
        best_tm = self.calculate_tm_nearest_neighbor(best_seq)
        
        # Try extending or shortening
        for length in range(self.min_annealing_length, self.max_annealing_length + 1):
            if direction == 'forward':
                test_end = min(start + length, len(template))
                test_seq = template[start:test_end]
            else:
                test_begin = max(0, start - length + 1)
                test_seq = self.reverse_complement(template[test_begin:start + 1])
            
            test_tm = self.calculate_tm_nearest_neighbor(test_seq)
            
            # Prefer sequences within target Tm range
            if self.target_tm_range[0] <= test_tm <= self.target_tm_range[1]:
                if not (self.target_tm_range[0] <= best_tm <= self.target_tm_range[1]):
                    best_seq = test_seq
                    best_tm = test_tm
                elif abs(test_tm - 60) < abs(best_tm - 60):  # Prefer closer to 60°C
                    best_seq = test_seq
                    best_tm = test_tm
        
        validation = self.validate_primer(best_seq)
        
        return best_seq, validation
    
    def design_gibson_primer(self, template: str, start: int, homology_seq: str,
                            direction: str = 'forward') -> Tuple[str, Dict]:
        """
        Design complete Gibson assembly primer with homology overhang + annealing region
        
        Args:
            template: DNA template sequence
            start: Starting position for annealing on template
            homology_seq: Homology sequence to add as 5' overhang
            direction: 'forward' or 'reverse'
        
        Returns:
            (complete_primer, validation_results)
        """
        # Design annealing region
        annealing_seq, anneal_validation = self.design_annealing_region(
            template, start, direction
        )
        
        # Combine homology + annealing
        if direction == 'forward':
            complete_primer = homology_seq.upper() + annealing_seq
        else:
            complete_primer = self.reverse_complement(homology_seq).upper() + annealing_seq
        
        # For Gibson primers, validate annealing region separately (not full primer)
        # Calculate annealing region metrics
        annealing_tm = self.calculate_tm_nearest_neighbor(annealing_seq)
        annealing_gc = self.calculate_gc_content(annealing_seq)
        
        # Calculate homology region Tm
        homology_seq_proper = homology_seq if direction == 'forward' else self.reverse_complement(homology_seq)
        homology_tm = self.calculate_tm_nearest_neighbor(homology_seq)
        
        # Validate only the annealing region
        annealing_validation = self.validate_primer(annealing_seq)
        
        # Create validation for full primer but with annealing-focused warnings
        full_validation = {
            'valid': annealing_validation['valid'],
            'warnings': [],
            'errors': annealing_validation['errors'].copy(),
            'metrics': {
                'length': len(complete_primer),
                'gc_content': self.calculate_gc_content(complete_primer),
                'tm': self.calculate_tm_nearest_neighbor(complete_primer),
                'has_gc_clamp': self.has_gc_clamp(annealing_seq),
                'hairpin_score': annealing_validation['metrics'].get('hairpin_score', 0),
                'self_dimer_score': annealing_validation['metrics'].get('self_dimer_score', 0)
            }
        }
        
        # Add region-specific info
        full_validation['annealing_region'] = {
            'sequence': annealing_seq,
            'length': len(annealing_seq),
            'tm': annealing_tm,
            'gc_content': annealing_gc,
            'validation': annealing_validation
        }
        full_validation['homology_region'] = {
            'sequence': homology_seq_proper,
            'length': len(homology_seq),
            'tm': homology_tm,
            'gc_content': self.calculate_gc_content(homology_seq)
        }
        
        # Add annealing region warnings (GC, Tm, GC clamp)
        if annealing_gc < self.gc_content_range[0] or annealing_gc > self.gc_content_range[1]:
            full_validation['warnings'].append(
                f"Annealing region GC content {annealing_gc}% is outside recommended range {self.gc_content_range[0]}-{self.gc_content_range[1]}%"
            )
        
        if annealing_tm < self.target_tm_range[0] or annealing_tm > self.target_tm_range[1]:
            full_validation['warnings'].append(
                f"Annealing region Tm {annealing_tm:.1f}°C is outside recommended range {self.target_tm_range[0]}-{self.target_tm_range[1]}°C"
            )
        
        if not full_validation['metrics']['has_gc_clamp']:
            full_validation['warnings'].append("Annealing region has no GC clamp at 3' end (recommend 1-2 G/C bases)")
        
        # Copy structure warnings from annealing validation (hairpin, dimers)
        for warning in annealing_validation['warnings']:
            if 'GC content' not in warning and 'Tm' not in warning and 'GC clamp' not in warning:
                full_validation['warnings'].append(warning)
        
        # Check homology region Tm (Takara recommends 50-60°C for Gibson assembly)
        if homology_tm < 50 or homology_tm > 60:
            full_validation['warnings'].append(
                f"Homology region Tm {homology_tm:.1f}°C is outside recommended range 50-60°C"
            )
        
        return complete_primer, full_validation
