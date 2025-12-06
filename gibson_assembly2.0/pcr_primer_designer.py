"""
PCR Primer Designer - Optimized primer design for standard PCR amplification
"""

from typing import Dict, List, Tuple
from primer_design import PrimerDesigner


class PCRPrimerOptimizer:
    """Optimize PCR primers for amplifying a target sequence"""
    
    def __init__(self):
        self.designer = PrimerDesigner()
    
    def design_pcr_primers(self, template_seq: str,
                          fwd_start: int = 0, fwd_end: int = None,
                          rev_start: int = None, rev_end: int = None,
                          min_primer_length: int = 18,
                          max_primer_length: int = 25,
                          target_tm: float = 60.0,
                          step: int = 1,
                          max_results: int = 100) -> List[Dict]:
        """
        Design PCR primers to amplify template sequence
        
        Args:
            template_seq: DNA template sequence
            fwd_start: Start position for forward primer range (default: 0)
            fwd_end: End position for forward primer range (default: 50)
            rev_start: Start position for reverse primer range (default: len-50)
            rev_end: End position for reverse primer range (default: len)
            min_primer_length: Minimum primer length
            max_primer_length: Maximum primer length
            target_tm: Target melting temperature
            step: Step size for scanning positions
            max_results: Maximum number of results to return
        
        Returns:
            List of primer pair results sorted by score
        """
        # Set defaults
        if fwd_end is None:
            fwd_end = min(50, len(template_seq) // 4)
        if rev_start is None:
            rev_start = max(len(template_seq) - 50, 3 * len(template_seq) // 4)
        if rev_end is None:
            rev_end = len(template_seq)
        
        # Validate ranges
        if fwd_start >= fwd_end:
            raise ValueError("Forward primer start must be less than end")
        if rev_start >= rev_end:
            raise ValueError("Reverse primer start must be less than end")
        if fwd_end >= rev_start:
            raise ValueError("Forward and reverse primer regions must not overlap")
        
        results = []
        total_tests = ((fwd_end - fwd_start) // step + 1) * ((rev_end - rev_start) // step + 1) * (max_primer_length - min_primer_length + 1) ** 2
        
        print(f"Scanning forward primers {fwd_start}-{fwd_end}, reverse primers {rev_start}-{rev_end}")
        print(f"Primer length range: {min_primer_length}-{max_primer_length} nt")
        print(f"Estimated tests: {total_tests:,}")
        
        if total_tests > 50000:
            print(f"⚠ WARNING: High complexity! Consider increasing step size or reducing ranges.")
        
        test_count = 0
        
        # Scan forward primer positions
        for fwd_pos in range(fwd_start, fwd_end + 1, step):
            # Scan reverse primer positions
            for rev_pos in range(rev_start, rev_end + 1, step):
                
                # Ensure amplicon is reasonable size
                amplicon_size = rev_pos - fwd_pos
                if amplicon_size < 50 or amplicon_size > 10000:
                    continue
                
                # Try different primer lengths
                for fwd_len in range(min_primer_length, max_primer_length + 1):
                    for rev_len in range(min_primer_length, max_primer_length + 1):
                        
                        try:
                            test_count += 1
                            if test_count % 1000 == 0:
                                print(f"  Progress: {test_count:,} tests...", end='\r')
                            
                            # Extract primers
                            if fwd_pos + fwd_len > len(template_seq):
                                continue
                            if rev_pos + rev_len > len(template_seq):
                                continue
                            
                            fwd_primer = template_seq[fwd_pos:fwd_pos + fwd_len]
                            rev_region = template_seq[rev_pos:rev_pos + rev_len]
                            
                            # Validate primers
                            fwd_val = self.designer.validate_primer(fwd_primer, is_gibson=False)
                            rev_val = self.designer.validate_primer(
                                self.designer.reverse_complement(rev_region), 
                                is_gibson=False
                            )
                            
                            # Calculate scores
                            fwd_score = self._score_primer(fwd_val, target_tm)
                            rev_score = self._score_primer(rev_val, target_tm)
                            
                            # Score pair compatibility
                            tm_diff = abs(fwd_val['metrics']['tm'] - rev_val['metrics']['tm'])
                            pair_score = 100
                            if tm_diff > 5:
                                pair_score -= 20 * (tm_diff - 5)
                            
                            overall_score = (fwd_score + rev_score + pair_score) / 3
                            
                            # Count issues
                            total_warnings = len(fwd_val['warnings']) + len(rev_val['warnings'])
                            total_errors = len(fwd_val['errors']) + len(rev_val['errors'])
                            
                            results.append({
                                'fwd_pos': fwd_pos,
                                'rev_pos': rev_pos,
                                'amplicon_size': amplicon_size,
                                'score': overall_score,
                                'fwd_score': fwd_score,
                                'rev_score': rev_score,
                                'pair_score': pair_score,
                                'primers': {
                                    'forward': fwd_primer,
                                    'reverse': self.designer.reverse_complement(rev_region)
                                },
                                'validations': {
                                    'forward': fwd_val,
                                    'reverse': rev_val
                                },
                                'warnings': total_warnings,
                                'errors': total_errors
                            })
                            
                            # Early termination
                            if max_results and len(results) >= max_results * 3:
                                results.sort(key=lambda x: x['score'], reverse=True)
                                results = results[:max_results * 2]
                        
                        except Exception as e:
                            continue
        
        print(f"\nCompleted {test_count:,} tests, found {len(results)} valid primer pairs")
        
        # Sort by score
        results.sort(key=lambda x: x['score'], reverse=True)
        
        return results[:max_results]
    
    def _score_primer(self, validation: Dict, target_tm: float) -> float:
        """Score a single primer"""
        score = 100
        metrics = validation['metrics']
        
        # Tm deviation from target
        tm_diff = abs(metrics['tm'] - target_tm)
        if tm_diff > 5:
            score -= 5 * (tm_diff - 5)
        
        # GC content
        gc = metrics['gc_content']
        if gc < 40 or gc > 60:
            score -= abs(50 - gc)
        
        # Penalties for warnings and errors
        score -= 10 * len(validation['warnings'])
        score -= 25 * len(validation['errors'])
        
        return max(0, min(100, score))
