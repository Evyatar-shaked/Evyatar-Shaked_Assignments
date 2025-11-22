"""
Example usage of the Gibson Assembly Primer Designer
"""

from primer_design import PrimerDesigner

def example_1_basic_validation():
    """Example 1: Basic primer validation"""
    print("=" * 70)
    print("EXAMPLE 1: Basic Primer Validation")
    print("=" * 70)
    
    designer = PrimerDesigner()
    
    # Test a simple primer
    primer = "GCTAGCTAGCTAGCTAGCTA"
    
    print(f"\nPrimer: 5'- {primer} -3'")
    
    # Calculate metrics
    tm = designer.calculate_tm_nearest_neighbor(primer)
    gc = designer.calculate_gc_content(primer)
    
    print(f"Length: {len(primer)} nt")
    print(f"Tm: {tm}°C")
    print(f"GC%: {gc}%")
    
    # Check for issues
    has_hairpin, score = designer.check_hairpin(primer)
    print(f"Hairpin: {'⚠️ Yes' if has_hairpin else '✅ No'} (score: {score})")


def example_2_gibson_primer_design():
    """Example 2: Design a complete Gibson primer"""
    print("\n\n" + "=" * 70)
    print("EXAMPLE 2: Gibson Primer Design")
    print("=" * 70)
    
    designer = PrimerDesigner()
    
    # Template sequence (e.g., start of a gene)
    template = "ATGGCTAGCATGACTGGTGGACAGCAAATGGGTCGCGGATCCGAATTCGAGCTCCGTCGAC"
    
    # Homology to add (e.g., from vector)
    homology = "CTGCAGGTCGACTCTAGAGG"
    
    print(f"\nTemplate: {template[:30]}...")
    print(f"Homology: {homology}")
    
    # Design forward primer
    primer, validation = designer.design_gibson_primer(
        template=template,
        start=0,
        homology_seq=homology,
        direction='forward'
    )
    
    print(f"\n✅ Gibson Primer Designed:")
    print(f"   5'- {primer} -3'")
    print(f"   Length: {len(primer)} nt")
    print(f"   Annealing Tm: {validation['annealing_region']['tm']}°C")
    print(f"   Homology Tm: {validation['homology_region']['tm']}°C")
    
    if validation['warnings']:
        print(f"\n⚠️  Warnings:")
        for warning in validation['warnings']:
            print(f"    • {warning}")


def example_3_primer_pair_compatibility():
    """Example 3: Check primer pair compatibility"""
    print("\n\n" + "=" * 70)
    print("EXAMPLE 3: Primer Pair Compatibility Check")
    print("=" * 70)
    
    designer = PrimerDesigner()
    
    # Design a primer pair
    template = "ATGGCTAGCATGACTGGTGGACAGCAAATGGGTCGCGGATCCGAATTCGAGCTCCGTCGAC"
    
    # Forward primer
    fwd_primer, _ = designer.design_annealing_region(
        template=template,
        start=0,
        direction='forward',
        target_length=20
    )
    
    # Reverse primer
    rev_primer, _ = designer.design_annealing_region(
        template=template,
        start=len(template) - 1,
        direction='reverse',
        target_length=20
    )
    
    print(f"\nForward: 5'- {fwd_primer} -3'")
    print(f"Reverse: 5'- {rev_primer} -3'")
    
    # Check Tm compatibility
    fwd_tm = designer.calculate_tm_nearest_neighbor(fwd_primer)
    rev_tm = designer.calculate_tm_nearest_neighbor(rev_primer)
    tm_diff = abs(fwd_tm - rev_tm)
    
    print(f"\nForward Tm: {fwd_tm}°C")
    print(f"Reverse Tm: {rev_tm}°C")
    print(f"Tm Difference: {tm_diff:.1f}°C", end="")
    
    if tm_diff <= 3:
        print(" ✅ (Good!)")
    else:
        print(" ⚠️  (Consider redesigning)")
    
    # Check for hetero-dimer
    has_hetero, score = designer.check_hetero_dimer(fwd_primer, rev_primer)
    print(f"Hetero-dimer: {'⚠️ Detected' if has_hetero else '✅ None'} (score: {score})")


def example_4_batch_validation():
    """Example 4: Validate multiple primers"""
    print("\n\n" + "=" * 70)
    print("EXAMPLE 4: Batch Primer Validation")
    print("=" * 70)
    
    designer = PrimerDesigner()
    
    primers = {
        "Primer_A": "GCTAGCTAGCTAGCTAGCTA",
        "Primer_B": "ATATATATATATATATATAT",  # Low GC, no clamp
        "Primer_C": "GCGCGCGCGCGCGCGCGCGC",  # High GC
        "Primer_D": "ACGTACGTACGTACGTACGT",  # Good
    }
    
    print("\nValidating primers...\n")
    
    for name, seq in primers.items():
        validation = designer.validate_primer(seq)
        
        status = "✅" if not validation['warnings'] else "⚠️"
        print(f"{status} {name}: Tm={validation['metrics']['tm']}°C, "
              f"GC={validation['metrics']['gc_content']}%, "
              f"Warnings={len(validation['warnings'])}")


def example_5_tm_calculation_methods():
    """Example 5: Compare Tm calculation methods"""
    print("\n\n" + "=" * 70)
    print("EXAMPLE 5: Tm Calculation Method Comparison")
    print("=" * 70)
    
    designer = PrimerDesigner()
    
    test_primers = [
        "ATCGATCGATCGATCG",
        "GCTAGCTAGCTAGCTA",
        "AAAATTTTAAAATTTT",
        "GCGCGCGCGCGCGCGC"
    ]
    
    print(f"\n{'Primer':<20} {'Wallace':<12} {'Nearest-Neighbor':<18} {'Difference'}")
    print("-" * 70)
    
    for primer in test_primers:
        tm_wallace = designer.calculate_tm_wallace(primer)
        tm_nn = designer.calculate_tm_nearest_neighbor(primer)
        diff = abs(tm_wallace - tm_nn)
        
        print(f"{primer:<20} {tm_wallace:<12.1f} {tm_nn:<18.1f} {diff:.1f}°C")
    
    print("\n💡 Note: Nearest-neighbor method is more accurate for long primers")


def main():
    """Run all examples"""
    print("\n" + "🧬" * 35)
    print("Gibson Assembly Primer Designer - Examples")
    print("🧬" * 35)
    
    example_1_basic_validation()
    example_2_gibson_primer_design()
    example_3_primer_pair_compatibility()
    example_4_batch_validation()
    example_5_tm_calculation_methods()
    
    print("\n\n" + "=" * 70)
    print("Examples completed! Try running gibson_wizard.py for interactive design.")
    print("=" * 70 + "\n")


if __name__ == "__main__":
    main()
