"""
Primer Optimization Module
Finds the best cut site in a vector for Gibson assembly by scoring primers
"""

from typing import Dict, List, Tuple, Optional
from primer_design import PrimerDesigner


class PrimerOptimizer:
    """Optimize Gibson assembly primer design by finding best cut sites"""
    
    def __init__(self):
        self.designer = PrimerDesigner()
    
    def score_primer(self, validation: Dict) -> float:
        """
        Score a primer based on validation results (higher is better)
        
        Args:
            validation: Validation dictionary from primer design
        
        Returns:
            Score (0-100, higher is better)
        """
        score = 100.0
        metrics = validation['metrics']
        
        # Tm scoring (prefer 55-65°C)
        if 'annealing_region' in validation:
            tm = validation['annealing_region']['tm']
        else:
            tm = metrics['tm']
        
        if tm < 50:
            score -= 30
        elif tm < 55:
            score -= 15
        elif tm > 65:
            score -= 10 * (tm - 65)
        
        # GC content scoring (prefer 40-60%)
        gc = metrics['gc_content']
        if gc < 30 or gc > 70:
            score -= 30
        elif gc < 40 or gc > 60:
            score -= 15
        
        # GC clamp bonus
        if metrics['has_gc_clamp']:
            score += 5
        else:
            score -= 5
        
        # Penalty for hairpins
        if metrics.get('hairpin_score', 0) >= 4:
            score -= 10 * metrics['hairpin_score']
        
        # Penalty for self-dimers
        if metrics.get('self_dimer_score', 0) >= 4:
            score -= 10 * metrics['self_dimer_score']
        
        # Penalty for hetero-dimers
        if metrics.get('hetero_dimer_score', 0) >= 4:
            score -= 10 * metrics['hetero_dimer_score']
        
        # Length penalty (prefer 40-60 nt total)
        length = metrics['length']
        if length < 35:
            score -= 5
        elif length > 70:
            score -= 10
        
        return max(0, min(100, score))
    
    def score_primer_pair(self, fwd_validation: Dict, rev_validation: Dict) -> float:
        """
        Score a primer pair (considers compatibility)
        
        Args:
            fwd_validation: Forward primer validation
            rev_validation: Reverse primer validation
        
        Returns:
            Combined score (0-100)
        """
        fwd_score = self.score_primer(fwd_validation)
        rev_score = self.score_primer(rev_validation)
        
        # Average of individual scores
        pair_score = (fwd_score + rev_score) / 2
        
        # Bonus/penalty for Tm compatibility
        if 'annealing_region' in fwd_validation and 'annealing_region' in rev_validation:
            fwd_tm = fwd_validation['annealing_region']['tm']
            rev_tm = rev_validation['annealing_region']['tm']
            tm_diff = abs(fwd_tm - rev_tm)
            
            if tm_diff <= 1:
                pair_score += 5
            elif tm_diff <= 3:
                pass  # No penalty
            elif tm_diff <= 5:
                pair_score -= 5
            else:
                pair_score -= 10 * (tm_diff - 5)
        
        return max(0, min(100, pair_score))
    
    def optimize_cut_site(self, vector_seq: str, insert_seq: str,
                         vector_upstream_start: int, vector_upstream_end: int,
                         vector_downstream_start: int, vector_downstream_end: int,
                         insert_fwd_start: int = None, insert_fwd_end: int = None,
                         insert_rev_start: int = None, insert_rev_end: int = None,
                         optimize_insert: bool = False,
                         homology_length: int = 25,
                         step: int = 10,
                         max_results: int = None) -> List[Dict]:
        """
        Find optimal primer positions by scoring all combinations
        
        Args:
            vector_seq: Vector DNA sequence
            insert_seq: Insert DNA sequence
            vector_upstream_start: Start of range for vector forward primer
            vector_upstream_end: End of range for vector forward primer
            vector_downstream_start: Start of range for vector reverse primer
            vector_downstream_end: End of range for vector reverse primer
            insert_fwd_start: Start of range for insert forward primer [optional]
            insert_fwd_end: End of range for insert forward primer [optional]
            insert_rev_start: Start of range for insert reverse primer [optional]
            insert_rev_end: End of range for insert reverse primer [optional]
            optimize_insert: Whether to optimize insert boundaries
            homology_length: Length of homology regions
            step: Step size for scanning (default 10 bp)
        
        Returns:
            List of results sorted by score (best first)
        """
        results = []
        
        # Ensure vector ranges are in correct order (but keep user values)
        # Swap if backwards
        if vector_upstream_start > vector_upstream_end:
            vector_upstream_start, vector_upstream_end = vector_upstream_end, vector_upstream_start
        
        # Swap if backwards
        if vector_downstream_start > vector_downstream_end:
            vector_downstream_start, vector_downstream_end = vector_downstream_end, vector_downstream_start
        
        # If insert ranges not specified or optimize_insert is False, use full insert
        if not optimize_insert:
            insert_fwd_start = insert_fwd_end = 0
            insert_rev_start = insert_rev_end = len(insert_seq)
        else:
            if insert_fwd_start is None:
                insert_fwd_start = 0
            if insert_fwd_end is None:
                insert_fwd_end = min(100, len(insert_seq) // 4)
            if insert_rev_start is None:
                insert_rev_start = max(len(insert_seq) - 100, 3 * len(insert_seq) // 4)
            if insert_rev_end is None:
                insert_rev_end = len(insert_seq)
            
            # Ensure insert ranges are in correct order (but keep user values)
            # Swap if backwards
            if insert_fwd_start > insert_fwd_end:
                insert_fwd_start, insert_fwd_end = insert_fwd_end, insert_fwd_start
            
            # Swap if backwards
            if insert_rev_start > insert_rev_end:
                insert_rev_start, insert_rev_end = insert_rev_end, insert_rev_start
        
        # Calculate total combinations
        vector_upstream_positions = ((vector_upstream_end - vector_upstream_start) // step + 1)
        vector_downstream_positions = ((vector_downstream_end - vector_downstream_start) // step + 1)
        total_tests = vector_upstream_positions * vector_downstream_positions
        
        if optimize_insert:
            insert_fwd_positions = ((insert_fwd_end - insert_fwd_start) // step + 1)
            insert_rev_positions = ((insert_rev_end - insert_rev_start) // step + 1)
            total_tests *= insert_fwd_positions * insert_rev_positions
            print(f"Scanning vector upstream {vector_upstream_start}-{vector_upstream_end}, downstream {vector_downstream_start}-{vector_downstream_end}")
            print(f"AND insert forward {insert_fwd_start}-{insert_fwd_end}, reverse {insert_rev_start}-{insert_rev_end} (step={step})...")
        else:
            print(f"Scanning vector upstream {vector_upstream_start}-{vector_upstream_end}, downstream {vector_downstream_start}-{vector_downstream_end} (step={step})...")
        
        print(f"Total combinations to test: {total_tests}")
        
        # Warn if complexity is very high
        if total_tests > 10000:
            print(f"⚠ WARNING: Testing {total_tests} combinations may take several minutes!")
            print(f"  Consider: increasing step size, reducing ranges, or disabling insert optimization")
        
        test_count = 0
        
        # Nested loops for all 4 ranges
        for vector_upstream_pos in range(vector_upstream_start, vector_upstream_end + 1, step):
            for vector_downstream_pos in range(vector_downstream_start, vector_downstream_end + 1, step):
                
                # Determine insert positions to test
                if optimize_insert:
                    insert_fwd_positions_range = range(insert_fwd_start, insert_fwd_end + 1, step)
                    insert_rev_positions_range = range(insert_rev_start, insert_rev_end + 1, step)
                else:
                    insert_fwd_positions_range = [0]
                    insert_rev_positions_range = [len(insert_seq)]
                
                for insert_fwd_pos in insert_fwd_positions_range:
                    for insert_rev_pos in insert_rev_positions_range:
                        
                        # Ensure insert region is valid
                        if insert_rev_pos - insert_fwd_pos < 20:
                            continue
                        
                        # Extract insert region to amplify
                        insert_region = insert_seq[insert_fwd_pos:insert_rev_pos]
                
                        try:
                            test_count += 1
                            if test_count % 100 == 0:
                                print(f"  Progress: {test_count}/{total_tests}...", end='\r')
                            
                            # Design vector primers
                            insert_start_homology = insert_region[:homology_length]
                            insert_end_homology = insert_region[-homology_length:]
                            
                            vector_fwd, vf_val = self.designer.design_gibson_primer(
                                template=vector_seq,
                                start=vector_upstream_pos,
                                homology_seq=insert_start_homology,
                                direction='forward'
                            )
                            
                            vector_rev, vr_val = self.designer.design_gibson_primer(
                                template=vector_seq,
                                start=vector_downstream_pos,
                                homology_seq=insert_end_homology,
                                direction='reverse'
                            )
                            
                            # Design insert primers
                            vector_upstream_homology = vector_seq[
                                max(0, vector_upstream_pos - homology_length):vector_upstream_pos
                            ]
                            vector_downstream_homology = vector_seq[
                                vector_downstream_pos:min(len(vector_seq), vector_downstream_pos + homology_length)
                            ]
                            
                            # Handle circular vector
                            if len(vector_upstream_homology) < homology_length:
                                needed = homology_length - len(vector_upstream_homology)
                                vector_upstream_homology = vector_seq[-needed:] + vector_upstream_homology
                            
                            if len(vector_downstream_homology) < homology_length:
                                needed = homology_length - len(vector_downstream_homology)
                                vector_downstream_homology = vector_downstream_homology + vector_seq[:needed]
                            
                            # Design insert primers for the specified region
                            insert_fwd, if_val = self.designer.design_gibson_primer(
                                template=insert_region,
                                start=0,
                                homology_seq=vector_upstream_homology,
                                direction='forward'
                            )
                            
                            insert_rev, ir_val = self.designer.design_gibson_primer(
                                template=insert_region,
                                start=len(insert_region) - 1,
                                homology_seq=vector_downstream_homology,
                                direction='reverse'
                            )
                            
                            # Score primers
                            vector_score = self.score_primer_pair(vf_val, vr_val)
                            insert_score = self.score_primer_pair(if_val, ir_val)
                            overall_score = (vector_score + insert_score) / 2
                            
                            # Count warnings and errors
                            total_warnings = sum([
                                len(vf_val['warnings']),
                                len(vr_val['warnings']),
                                len(if_val['warnings']),
                                len(ir_val['warnings'])
                            ])
                            
                            total_errors = sum([
                                len(vf_val['errors']),
                                len(vr_val['errors']),
                                len(if_val['errors']),
                                len(ir_val['errors'])
                            ])
                            
                            # Generate final construct
                            from primer_design import generate_final_construct
                            final_construct = generate_final_construct(
                                vector_seq, insert_seq, vector_upstream_pos, insert_fwd_pos, insert_rev_pos
                            )
                            
                            results.append({
                                'vector_upstream_pos': vector_upstream_pos,
                                'vector_downstream_pos': vector_downstream_pos,
                                'insert_fwd_pos': insert_fwd_pos,
                                'insert_rev_pos': insert_rev_pos,
                                'score': overall_score,
                                'vector_score': vector_score,
                                'insert_score': insert_score,
                                'primers': {
                                    'vector_forward': vector_fwd,
                                    'vector_reverse': vector_rev,
                                    'insert_forward': insert_fwd,
                                    'insert_reverse': insert_rev
                                },
                                'validations': {
                                    'vector_forward': vf_val,
                                    'vector_reverse': vr_val,
                                    'insert_forward': if_val,
                                    'insert_reverse': ir_val
                                },
                                'warnings': total_warnings,
                                'errors': total_errors,
                                'final_construct': final_construct
                            })
                            
                            # Early termination if we have enough high-quality results
                            if max_results and len(results) >= max_results * 3:
                                # Keep only best results so far
                                results.sort(key=lambda x: x['score'], reverse=True)
                                results = results[:max_results * 2]
                        
                        except Exception as e:
                            # Skip positions that cause errors
                            continue
        
        # Sort by score (highest first)
        results.sort(key=lambda x: x['score'], reverse=True)
        
        return results
    
    def find_best_cut_site(self, vector_seq: str, insert_seq: str,
                          start_range: int, end_range: int,
                          homology_length: int = 25) -> Dict:
        """
        Find the single best cut site (convenience method)
        
        Returns:
            Best result dictionary
        """
        results = self.optimize_cut_site(
            vector_seq, insert_seq, start_range, end_range, homology_length
        )
        
        if results:
            return results[0]
        else:
            raise ValueError("No valid cut sites found in range")
    
    def get_top_n_results(self, vector_seq: str, insert_seq: str,
                         vector_upstream_start: int = None, vector_upstream_end: int = None,
                         vector_downstream_start: int = None, vector_downstream_end: int = None,
                         insert_fwd_start: int = None, insert_fwd_end: int = None,
                         insert_rev_start: int = None, insert_rev_end: int = None,
                         optimize_insert: bool = False,
                         n: int = 5, homology_length: int = 25, step: int = 10) -> List[Dict]:
        """
        Get top N best primer combinations
        
        Args:
            vector_seq: Vector DNA sequence
            insert_seq: Insert DNA sequence
            vector_upstream_start: Start of range for vector forward primer
            vector_upstream_end: End of range for vector forward primer
            vector_downstream_start: Start of range for vector reverse primer
            vector_downstream_end: End of range for vector reverse primer
            insert_fwd_start: Start of range for insert forward primer (optional)
            insert_fwd_end: End of range for insert forward primer (optional)
            insert_rev_start: Start of range for insert reverse primer (optional)
            insert_rev_end: End of range for insert reverse primer (optional)
            optimize_insert: Whether to optimize insert boundaries
            n: Number of top results to return
            homology_length: Length of homology regions
            step: Step size for scanning
        
        Returns:
            List of top N results
        """
        # Set defaults if not provided
        if vector_upstream_start is None:
            vector_upstream_start = 100
        if vector_upstream_end is None:
            vector_upstream_end = min(500, len(vector_seq) // 3)
        if vector_downstream_start is None:
            vector_downstream_start = max(vector_upstream_end + 100, 2 * len(vector_seq) // 3)
        if vector_downstream_end is None:
            vector_downstream_end = len(vector_seq) - 100
        
        results = self.optimize_cut_site(
            vector_seq, insert_seq, 
            vector_upstream_start, vector_upstream_end,
            vector_downstream_start, vector_downstream_end,
            insert_fwd_start, insert_fwd_end,
            insert_rev_start, insert_rev_end,
            optimize_insert, homology_length, step,
            max_results=n
        )
        
        return results[:n]
