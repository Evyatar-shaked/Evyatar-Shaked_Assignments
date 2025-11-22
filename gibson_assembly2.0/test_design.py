"""
Simple test script for Gibson Assembly Primer Designer
Demonstrates a complete binary assembly design without user interaction
"""

from primer_design import PrimerDesigner

def main():
    print("=" * 80)
    print("AUTOMATED GIBSON ASSEMBLY PRIMER DESIGN TEST")
    print("=" * 80)
    
    # Initialize designer
    designer = PrimerDesigner()
    
    # Test sequences
    # Small vector (pUC19-like backbone, 500 bp excerpt)
    vector = """TTCTCATGTTTGACAGCTTATCATCGATAAGCTTTAATGCGGTAGTTTATCACAGTTAAATTGCTAACGCAGTCAGG
CACCGTGTATGAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTCACCCTGGATGCTGTAGGCATAGGCTTGG
TTATGCCGGTACTGCCGGGCCTCTTGCGGGATATCGTCCATTCCGACAGCATCGCCAGTCACTATGGCGTGCTGCTAG
CGCTATATGCGTTGATGCAATTTCTATGCGCACCCGTTCTCGGAGCACTGTCCGACCGCTTTGGCCGCCGCCCAGTCC
TGCTCGCTTCGCTACTTGGAGCCACTATCGACTACGCGATCATGGCGACCACACCCGTCCTGTGGATCCTCTACGCCG
GACGCATCGTGGCCGGCATCACCGGCGCCACAGGTGCGGTTGCTGGCGCCTATATCGCCGACATCACCGATGGGGAAG
ATCGGGCTCGCCACTTCGGGCTCATGAGCGCTTGTTTCGGCGTGGGTATGGTGGCAGGCCCCGTGGCC""".replace("\n", "")
    
    # GFP insert (300 bp excerpt)
    insert = """ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCC
ACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCG
GCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACC
ACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACG
ACGGCAACTACAAGACCCGCGCCGAG""".replace("\n", "")
    
    vector_cut_site = 250  # Middle of vector
    homology_length = 25
    
    print(f"\nTest Configuration:")
    print(f"  Vector length: {len(vector)} bp")
    print(f"  Insert length: {len(insert)} bp")
    print(f"  Cut site: position {vector_cut_site + 1}")
    print(f"  Homology length: {homology_length} nt")
    
    print("\n" + "-" * 80)
    print("DESIGNING PRIMERS...")
    print("-" * 80)
    
    # Design vector primers
    print("\n1. Vector Forward Primer (linearize vector, add insert start homology)")
    insert_start_homology = insert[:homology_length]
    
    vector_fwd, vf_val = designer.design_gibson_primer(
        template=vector,
        start=max(0, vector_cut_site - 50),
        homology_seq=insert_start_homology,
        direction='forward'
    )
    
    print(f"   Sequence: 5'- {vector_fwd} -3'")
    print(f"   Length: {len(vector_fwd)} nt")
    print(f"   Annealing Tm: {vf_val['annealing_region']['tm']}°C")
    print(f"   Status: {'✅ PASS' if not vf_val['errors'] else '❌ FAIL'}")
    
    print("\n2. Vector Reverse Primer (linearize vector, add insert end homology)")
    insert_end_homology = insert[-homology_length:]
    
    vector_rev, vr_val = designer.design_gibson_primer(
        template=vector,
        start=min(len(vector) - 1, vector_cut_site + 50),
        homology_seq=insert_end_homology,
        direction='reverse'
    )
    
    print(f"   Sequence: 5'- {vector_rev} -3'")
    print(f"   Length: {len(vector_rev)} nt")
    print(f"   Annealing Tm: {vr_val['annealing_region']['tm']}°C")
    print(f"   Status: {'✅ PASS' if not vr_val['errors'] else '❌ FAIL'}")
    
    # Check vector primer compatibility
    vf_tm = vf_val['annealing_region']['tm']
    vr_tm = vr_val['annealing_region']['tm']
    v_tm_diff = abs(vf_tm - vr_tm)
    print(f"\n   Vector primer Tm difference: {v_tm_diff:.1f}°C", end="")
    if v_tm_diff <= 3:
        print(" ✅")
    else:
        print(" ⚠️")
    
    # Design insert primers
    print("\n3. Insert Forward Primer (amplify insert, add vector upstream homology)")
    vector_upstream_homology = vector[max(0, vector_cut_site - homology_length):vector_cut_site]
    if len(vector_upstream_homology) < homology_length:
        needed = homology_length - len(vector_upstream_homology)
        vector_upstream_homology = vector[-needed:] + vector_upstream_homology
    
    insert_fwd, if_val = designer.design_gibson_primer(
        template=insert,
        start=0,
        homology_seq=vector_upstream_homology,
        direction='forward'
    )
    
    print(f"   Sequence: 5'- {insert_fwd} -3'")
    print(f"   Length: {len(insert_fwd)} nt")
    print(f"   Annealing Tm: {if_val['annealing_region']['tm']}°C")
    print(f"   Status: {'✅ PASS' if not if_val['errors'] else '❌ FAIL'}")
    
    print("\n4. Insert Reverse Primer (amplify insert, add vector downstream homology)")
    vector_downstream_homology = vector[vector_cut_site:min(len(vector), vector_cut_site + homology_length)]
    if len(vector_downstream_homology) < homology_length:
        needed = homology_length - len(vector_downstream_homology)
        vector_downstream_homology = vector_downstream_homology + vector[:needed]
    
    insert_rev, ir_val = designer.design_gibson_primer(
        template=insert,
        start=len(insert) - 1,
        homology_seq=vector_downstream_homology,
        direction='reverse'
    )
    
    print(f"   Sequence: 5'- {insert_rev} -3'")
    print(f"   Length: {len(insert_rev)} nt")
    print(f"   Annealing Tm: {ir_val['annealing_region']['tm']}°C")
    print(f"   Status: {'✅ PASS' if not ir_val['errors'] else '❌ FAIL'}")
    
    # Check insert primer compatibility
    if_tm = if_val['annealing_region']['tm']
    ir_tm = ir_val['annealing_region']['tm']
    i_tm_diff = abs(if_tm - ir_tm)
    print(f"\n   Insert primer Tm difference: {i_tm_diff:.1f}°C", end="")
    if i_tm_diff <= 3:
        print(" ✅")
    else:
        print(" ⚠️")
    
    # Summary
    print("\n" + "=" * 80)
    print("VALIDATION SUMMARY")
    print("=" * 80)
    
    all_pass = True
    total_warnings = 0
    
    for name, val in [("Vector Forward", vf_val), ("Vector Reverse", vr_val),
                      ("Insert Forward", if_val), ("Insert Reverse", ir_val)]:
        warnings = len(val['warnings'])
        errors = len(val['errors'])
        total_warnings += warnings
        
        status = "✅ PASS" if errors == 0 else "❌ FAIL"
        print(f"{name:20} {status:10} Warnings: {warnings}, Errors: {errors}")
        
        if errors > 0:
            all_pass = False
    
    print("\n" + "=" * 80)
    if all_pass:
        print("✅ ALL PRIMERS DESIGNED SUCCESSFULLY!")
        if total_warnings > 0:
            print(f"⚠️  {total_warnings} warning(s) detected - review before ordering")
        else:
            print("🎉 NO WARNINGS - Primers are ready to order!")
    else:
        print("❌ SOME PRIMERS FAILED VALIDATION - Review and redesign")
    print("=" * 80)
    
    # PCR suggestions
    print("\n📋 SUGGESTED PCR CONDITIONS")
    print("-" * 80)
    print(f"Vector PCR:")
    print(f"  Annealing temperature: {(vf_tm + vr_tm) / 2 - 5:.1f}°C")
    print(f"  Expected product: ~{len(vector)} bp")
    
    print(f"\nInsert PCR:")
    print(f"  Annealing temperature: {(if_tm + ir_tm) / 2 - 5:.1f}°C")
    print(f"  Expected product: ~{len(insert)} bp")
    
    print("\n✨ Test completed successfully!")


if __name__ == "__main__":
    main()
