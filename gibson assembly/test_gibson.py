"""
Test script for gibson_primer_design module
Run this to verify the installation and basic functionality
"""

try:
    from gibson_primer_design import GibsonPrimerDesigner
    from pydna.dseqrecord import Dseqrecord
    print("✓ Imports successful!")
except ImportError as e:
    print(f"✗ Import error: {e}")
    print("\nPlease install required packages:")
    print("  pip install pydna biopython")
    exit(1)


def test_basic_design():
    """Test basic primer design functionality"""
    print("\n" + "="*80)
    print("TEST 1: Basic Primer Design")
    print("="*80)
    
    # Create simple test sequences
    vector_seq = "ATGCGTACGTTAGCCTAGGCTAGCTAGGCTTACG" * 20
    vector = Dseqrecord(vector_seq, circular=True, name="pTestVector")
    
    insert_seq = "GGGAAAATTTCCCGGGTTTAAACCCGGGAAAATTTCCCGGGTTTAAA"
    insert = Dseqrecord(insert_seq, name="TestInsert")
    
    print(f"Vector: {len(vector)} bp (circular)")
    print(f"Insert: {len(insert)} bp")
    print(f"Insert site: 100")
    
    # Design primers
    designer = GibsonPrimerDesigner(homology_length=40, target_tm=60.0)
    result = designer.design_primers(vector, insert, insert_site=100)
    
    # Check results
    assert 'primers' in result, "Missing primers in result"
    assert 'insert_fwd' in result['primers'], "Missing insert_fwd primer"
    assert 'insert_rev' in result['primers'], "Missing insert_rev primer"
    assert 'vector_fwd' in result['primers'], "Missing vector_fwd primer"
    assert 'vector_rev' in result['primers'], "Missing vector_rev primer"
    
    print("\n✓ All 4 primers generated successfully!")
    
    # Check primer lengths
    for name, seq in result['primers'].items():
        print(f"  {name}: {len(seq)} bp")
    
    # Check assembly
    if result['assembly_result']:
        print(f"\n✓ Assembly verification successful!")
        print(f"  Final plasmid: {len(result['assembly_result'])} bp")
    else:
        print("\n✗ Assembly verification failed")
        return False
    
    return True


def test_multi_step():
    """Test multi-step assembly planning"""
    print("\n" + "="*80)
    print("TEST 2: Multi-Step Assembly Planning")
    print("="*80)
    
    # Create test sequences
    vector = Dseqrecord("ATGCGTACGTTAGCCTAGGC" * 30, circular=True, name="pVector")
    inserts = [
        Dseqrecord("GGGAAAATTTCCCGGGTTTAAA", name="Insert1"),
        Dseqrecord("TTTGGGCCCAAAATTTGGGCCC", name="Insert2"),
    ]
    
    print(f"Vector: {len(vector)} bp")
    print(f"Inserts: {len(inserts)} sequences")
    
    designer = GibsonPrimerDesigner(homology_length=30)
    steps = designer.plan_multi_step_assembly(vector, inserts, insert_site=50)
    
    assert len(steps) == len(inserts), "Wrong number of steps"
    
    print(f"\n✓ Planned {len(steps)} assembly steps successfully!")
    
    for i, step in enumerate(steps, 1):
        print(f"\n  Step {i}:")
        print(f"    Input vector: {len(step['input_vector'])} bp")
        if step['assembly_result']:
            print(f"    Output vector: {len(step['assembly_result'])} bp")
        print(f"    Primers: {len(step['primers'])} designed")
    
    return True


def test_primer_details():
    """Test detailed primer information"""
    print("\n" + "="*80)
    print("TEST 3: Primer Details and Tm Calculation")
    print("="*80)
    
    vector = Dseqrecord("ATGCGTACGTTAGCCTAGGC" * 25, circular=True, name="pVector")
    insert = Dseqrecord("GGGAAAATTTCCCGGGTTTAAACCCGGG", name="Insert")
    
    designer = GibsonPrimerDesigner(homology_length=35, target_tm=60.0, tm_tolerance=5.0)
    result = designer.design_primers(vector, insert, insert_site=75)
    
    print("\nChecking primer details...")
    
    for primer_name in ['insert_fwd', 'insert_rev', 'vector_fwd', 'vector_rev']:
        details = result['primer_details'][primer_name]
        
        assert 'sequence' in details, f"Missing sequence for {primer_name}"
        assert 'length' in details, f"Missing length for {primer_name}"
        assert 'anneal_tm' in details, f"Missing Tm for {primer_name}"
        assert 'homology_length' in details, f"Missing homology_length for {primer_name}"
        
        print(f"\n  {primer_name}:")
        print(f"    Length: {details['length']} bp")
        print(f"    Homology: {details['homology_length']} bp")
        print(f"    Annealing Tm: {details['anneal_tm']}°C")
        print(f"    Full Tm: {details['full_tm']}°C")
    
    print("\n✓ All primer details validated!")
    return True


def test_circular_handling():
    """Test circular DNA handling at wrap-around point"""
    print("\n" + "="*80)
    print("TEST 4: Circular DNA Wrap-Around")
    print("="*80)
    
    # Small circular vector, insert near end to test wrap-around
    vector = Dseqrecord("ATGCGTACGTTAGCCTAGGC" * 10, circular=True, name="pSmall")
    insert = Dseqrecord("GGGAAAATTTCCCGGG", name="Insert")
    
    # Insert near the end (will wrap around for homology extraction)
    insert_site = len(vector) - 10
    
    print(f"Vector: {len(vector)} bp (circular)")
    print(f"Insert site: {insert_site} (near end)")
    print(f"This tests circular wrap-around for homology regions")
    
    designer = GibsonPrimerDesigner(homology_length=30)
    result = designer.design_primers(vector, insert, insert_site=insert_site)
    
    # Should still generate primers even at wrap-around point
    assert result['primers']['insert_fwd'], "Failed at wrap-around"
    print("\n✓ Circular wrap-around handled correctly!")
    
    return True


def main():
    """Run all tests"""
    print("\n" + "#"*80)
    print("# Gibson Primer Design Module - Test Suite")
    print("#"*80)
    
    tests = [
        ("Basic Primer Design", test_basic_design),
        ("Multi-Step Assembly", test_multi_step),
        ("Primer Details", test_primer_details),
        ("Circular DNA Handling", test_circular_handling),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"\n✗ Test failed with error: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))
    
    # Summary
    print("\n" + "#"*80)
    print("# TEST SUMMARY")
    print("#"*80)
    
    for name, success in results:
        status = "✓ PASSED" if success else "✗ FAILED"
        print(f"{status}: {name}")
    
    all_passed = all(success for _, success in results)
    
    if all_passed:
        print("\n" + "="*80)
        print("🎉 ALL TESTS PASSED! Module is working correctly.")
        print("="*80)
        print("\nYou can now use the module in your projects!")
        print("See gibson_demo.ipynb for detailed examples.")
    else:
        print("\n" + "="*80)
        print("⚠️  SOME TESTS FAILED")
        print("="*80)
        return False
    
    return True


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
