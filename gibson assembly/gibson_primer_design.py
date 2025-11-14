"""
Gibson Assembly Primer Design Module for Binary Cloning

This module uses pydna to automatically design primers for binary Gibson assembly,
where each step adds one insert to a vector. For each assembly step, 4 primers are designed:
- 2 primers for amplifying the insert (with vector homology tails)
- 2 primers for amplifying the vector (with insert homology tails)
"""

from typing import List, Dict, Tuple, Optional
from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pydna.amplify import pcr
from pydna.assembly import Assembly
from Bio.SeqUtils import MeltingTemp as mt


class GibsonPrimerDesigner:
    """
    Design primers for binary Gibson assembly using pydna.
    
    Each step adds one insert to a vector backbone, requiring 4 primers:
    1. Insert Forward: amplifies insert with vector left homology tail
    2. Insert Reverse: amplifies insert with vector right homology tail
    3. Vector Forward: amplifies vector with insert left homology tail
    4. Vector Reverse: amplifies vector with insert right homology tail
    """
    
    def __init__(self, homology_length: int = 40, min_anneal_length: int = 20, 
                 target_tm: float = 60.0, tm_tolerance: float = 5.0):
        """
        Initialize the primer designer.
        
        Args:
            homology_length: Length of homology region for Gibson assembly (default 40 bp)
            min_anneal_length: Minimum annealing length for primers (default 20 bp)
            target_tm: Target melting temperature for annealing region (default 60°C)
            tm_tolerance: Acceptable Tm deviation (default ±5°C)
        """
        self.homology_length = homology_length
        self.min_anneal_length = min_anneal_length
        self.target_tm = target_tm
        self.tm_tolerance = tm_tolerance
    
    def _calculate_tm(self, seq: str) -> float:
        """Calculate melting temperature using Tm_NN method."""
        return mt.Tm_NN(seq, Na=50, dnac1=250, dnac2=250)
    
    def _find_optimal_anneal_length(self, template_seq: str, start_pos: int, 
                                   direction: str = 'forward') -> Tuple[str, int, float]:
        """
        Find optimal annealing length based on Tm.
        
        Args:
            template_seq: DNA sequence to anneal to
            start_pos: Starting position for annealing
            direction: 'forward' or 'reverse' for annealing direction
            
        Returns:
            Tuple of (annealing_sequence, length, Tm)
        """
        max_length = min(40, len(template_seq))  # Don't go beyond 40 bp
        
        for length in range(self.min_anneal_length, max_length + 1):
            if direction == 'forward':
                anneal_seq = template_seq[start_pos:start_pos + length]
            else:
                anneal_seq = template_seq[max(0, start_pos - length):start_pos]
            
            if len(anneal_seq) < self.min_anneal_length:
                continue
                
            tm = self._calculate_tm(anneal_seq)
            
            if abs(tm - self.target_tm) <= self.tm_tolerance:
                return anneal_seq, length, tm
        
        # If no optimal found, use min_anneal_length
        if direction == 'forward':
            anneal_seq = template_seq[start_pos:start_pos + self.min_anneal_length]
        else:
            anneal_seq = template_seq[max(0, start_pos - self.min_anneal_length):start_pos]
        
        tm = self._calculate_tm(anneal_seq)
        return anneal_seq, len(anneal_seq), tm
    
    def design_primers(self, vector: Dseqrecord, insert: Dseqrecord, 
                      insert_site: int) -> Dict:
        """
        Design 4 primers for binary Gibson assembly.
        
        Args:
            vector: Vector backbone (pydna Dseqrecord, can be circular)
            insert: Insert sequence (pydna Dseqrecord)
            insert_site: Position in vector where insert will be placed (0-based)
            
        Returns:
            Dictionary containing:
                - primers: Dict with 'insert_fwd', 'insert_rev', 'vector_fwd', 'vector_rev'
                - primer_objects: pydna Primer objects
                - pcr_products: Simulated PCR products
                - assembly_result: Predicted assembly result
                - primer_details: Detailed info about each primer (Tm, lengths, etc.)
        """
        # Ensure vector is circular for proper handling
        if not vector.circular:
            print("Warning: Vector is not marked as circular. Setting circular=True.")
            vector = vector.looped()
        
        # Get vector sequence as string
        vec_seq = str(vector.seq)
        ins_seq = str(insert.seq)
        
        # Extract homology regions from vector (circular-aware)
        vec_left_homology = self._get_circular_subseq(vec_seq, 
                                                       insert_site - self.homology_length, 
                                                       self.homology_length)
        vec_right_homology = self._get_circular_subseq(vec_seq, 
                                                        insert_site, 
                                                        self.homology_length)
        
        # Design insert primers with optimal annealing length
        ins_fwd_anneal, ins_fwd_anneal_len, ins_fwd_tm = self._find_optimal_anneal_length(
            ins_seq, 0, 'forward'
        )
        ins_rev_anneal, ins_rev_anneal_len, ins_rev_tm = self._find_optimal_anneal_length(
            ins_seq, len(ins_seq), 'reverse'
        )
        
        # Insert forward: vector left homology + insert forward annealing
        insert_fwd_seq = vec_left_homology + ins_fwd_anneal
        
        # Insert reverse: RC(vector right homology) + RC(insert reverse annealing)
        insert_rev_seq = str(Dseqrecord(vec_right_homology).reverse_complement().seq) + \
                        str(Dseqrecord(ins_rev_anneal).reverse_complement().seq)
        
        # Design vector primers with optimal annealing length
        vec_fwd_anneal, vec_fwd_anneal_len, vec_fwd_tm = self._find_optimal_anneal_length(
            vec_seq, insert_site, 'reverse'  # Goes backwards from cut site
        )
        vec_rev_anneal, vec_rev_anneal_len, vec_rev_tm = self._find_optimal_anneal_length(
            vec_seq, insert_site, 'forward'  # Goes forward from cut site
        )
        
        # Vector forward: RC(insert left) + vector left annealing
        ins_left_homology = ins_seq[:self.homology_length]
        vector_fwd_seq = str(Dseqrecord(ins_left_homology).reverse_complement().seq) + vec_fwd_anneal
        
        # Vector reverse: insert right + RC(vector right annealing)
        ins_right_homology = ins_seq[-self.homology_length:]
        vector_rev_seq = ins_right_homology + \
                        str(Dseqrecord(vec_rev_anneal).reverse_complement().seq)
        
        # Create Primer objects
        insert_fwd_primer = Primer(insert_fwd_seq, name="Insert_Fwd")
        insert_rev_primer = Primer(insert_rev_seq, name="Insert_Rev")
        vector_fwd_primer = Primer(vector_fwd_seq, name="Vector_Fwd")
        vector_rev_primer = Primer(vector_rev_seq, name="Vector_Rev")
        
        # Simulate PCR products
        insert_pcr = pcr(insert_fwd_primer, insert_rev_primer, insert)
        vector_pcr = pcr(vector_fwd_primer, vector_rev_primer, vector)
        
        # Simulate Gibson assembly
        assembly = Assembly([insert_pcr, vector_pcr], limit=self.homology_length - 5)
        assembly_products = assembly.assemble_circular()
        
        # Prepare detailed results
        result = {
            'primers': {
                'insert_fwd': insert_fwd_seq,
                'insert_rev': insert_rev_seq,
                'vector_fwd': vector_fwd_seq,
                'vector_rev': vector_rev_seq,
            },
            'primer_objects': {
                'insert_fwd': insert_fwd_primer,
                'insert_rev': insert_rev_primer,
                'vector_fwd': vector_fwd_primer,
                'vector_rev': vector_rev_primer,
            },
            'pcr_products': {
                'insert': insert_pcr,
                'vector': vector_pcr,
            },
            'assembly_result': assembly_products[0] if assembly_products else None,
            'primer_details': {
                'insert_fwd': {
                    'sequence': insert_fwd_seq,
                    'length': len(insert_fwd_seq),
                    'homology_tail': vec_left_homology,
                    'homology_length': len(vec_left_homology),
                    'anneal_seq': ins_fwd_anneal,
                    'anneal_length': ins_fwd_anneal_len,
                    'anneal_tm': round(ins_fwd_tm, 1),
                    'full_tm': round(self._calculate_tm(insert_fwd_seq), 1),
                },
                'insert_rev': {
                    'sequence': insert_rev_seq,
                    'length': len(insert_rev_seq),
                    'homology_tail': str(Dseqrecord(vec_right_homology).reverse_complement().seq),
                    'homology_length': len(vec_right_homology),
                    'anneal_seq': str(Dseqrecord(ins_rev_anneal).reverse_complement().seq),
                    'anneal_length': ins_rev_anneal_len,
                    'anneal_tm': round(ins_rev_tm, 1),
                    'full_tm': round(self._calculate_tm(insert_rev_seq), 1),
                },
                'vector_fwd': {
                    'sequence': vector_fwd_seq,
                    'length': len(vector_fwd_seq),
                    'homology_tail': str(Dseqrecord(ins_left_homology).reverse_complement().seq),
                    'homology_length': len(ins_left_homology),
                    'anneal_seq': vec_fwd_anneal,
                    'anneal_length': vec_fwd_anneal_len,
                    'anneal_tm': round(vec_fwd_tm, 1),
                    'full_tm': round(self._calculate_tm(vector_fwd_seq), 1),
                },
                'vector_rev': {
                    'sequence': vector_rev_seq,
                    'length': len(vector_rev_seq),
                    'homology_tail': ins_right_homology,
                    'homology_length': len(ins_right_homology),
                    'anneal_seq': str(Dseqrecord(vec_rev_anneal).reverse_complement().seq),
                    'anneal_length': vec_rev_anneal_len,
                    'anneal_tm': round(vec_rev_tm, 1),
                    'full_tm': round(self._calculate_tm(vector_rev_seq), 1),
                },
            }
        }
        
        return result
    
    def _get_circular_subseq(self, seq: str, start: int, length: int) -> str:
        """Get subsequence handling circular wrapping."""
        n = len(seq)
        if length <= 0:
            return ''
        start = start % n
        if start + length <= n:
            return seq[start:start + length]
        else:
            return seq[start:] + seq[:(start + length) % n]
    
    def plan_multi_step_assembly(self, initial_vector: Dseqrecord, 
                                inserts: List[Dseqrecord], 
                                insert_site: int) -> List[Dict]:
        """
        Plan multi-step binary Gibson assembly.
        
        Args:
            initial_vector: Starting vector backbone
            inserts: List of insert sequences to add sequentially
            insert_site: Position where each insert will be added
            
        Returns:
            List of dictionaries, one per assembly step, each containing
            primer designs and assembly predictions
        """
        steps = []
        current_vector = initial_vector
        
        for i, insert in enumerate(inserts, start=1):
            print(f"\n{'='*60}")
            print(f"Planning Step {i}: Adding insert of length {len(insert)} bp")
            print(f"{'='*60}")
            
            result = self.design_primers(current_vector, insert, insert_site)
            result['step'] = i
            result['insert'] = insert
            result['input_vector'] = current_vector
            
            steps.append(result)
            
            # Update vector for next step
            if result['assembly_result']:
                current_vector = result['assembly_result']
                print(f"✓ Assembly successful. New vector length: {len(current_vector)} bp")
            else:
                print(f"✗ Warning: Assembly prediction failed for step {i}")
                break
        
        return steps
    
    def print_primer_summary(self, result: Dict):
        """Print a nicely formatted summary of the primer design."""
        print("\n" + "="*80)
        print("GIBSON ASSEMBLY PRIMER DESIGN SUMMARY")
        print("="*80)
        
        for primer_name in ['insert_fwd', 'insert_rev', 'vector_fwd', 'vector_rev']:
            details = result['primer_details'][primer_name]
            print(f"\n{primer_name.upper().replace('_', ' ')}:")
            print(f"  Sequence (5'→3'): {details['sequence']}")
            print(f"  Total Length: {details['length']} bp")
            print(f"  Homology Tail: {details['homology_length']} bp (for Gibson assembly)")
            print(f"  Annealing Region: {details['anneal_length']} bp (Tm: {details['anneal_tm']}°C)")
            print(f"  Full Primer Tm: {details['full_tm']}°C")
        
        print("\n" + "="*80)
        print("PCR PRODUCTS:")
        print("="*80)
        print(f"Insert PCR product: {len(result['pcr_products']['insert'])} bp")
        print(f"Vector PCR product: {len(result['pcr_products']['vector'])} bp")
        
        if result['assembly_result']:
            print(f"\n✓ Assembly Product: {len(result['assembly_result'])} bp")
        else:
            print(f"\n✗ Assembly failed - check homology regions")


def example_usage():
    """Example demonstrating the Gibson primer designer."""
    
    # Create example sequences
    # Vector: Simple circular plasmid
    vector_seq = "ATGCGTACGTTAGCCTAGGCTAGCTAGGCTTACGATGCGTACGTTAGCCTAGGCTAGCTAGGCTTACG" * 10
    vector = Dseqrecord(vector_seq, circular=True, name="pVector")
    
    # Insert: Gene to be cloned
    insert_seq = "GGGAAAATTTCCCGGGTTTAAACCCGGGAAAATTTCCCGGGTTTAAA"
    insert = Dseqrecord(insert_seq, name="Insert1")
    
    # Design primers
    designer = GibsonPrimerDesigner(homology_length=40, target_tm=60.0)
    result = designer.design_primers(vector, insert, insert_site=100)
    
    # Print summary
    designer.print_primer_summary(result)
    
    return result


if __name__ == "__main__":
    result = example_usage()
