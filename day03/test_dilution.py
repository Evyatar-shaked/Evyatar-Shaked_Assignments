"""
Comprehensive test suite for dilution calculator.
Tests the business logic, unit conversions, and edge cases.

Run with: python -m pytest test_dilution.py -v
Or: python test_dilution.py
"""

import unittest
import sys
from dilution_core import (
    calculate_dilution,
    convert_volume,
    convert_concentration,
    validate_dilution,
    VOLUME_UNITS,
    CONCENTRATION_UNITS,
    MOLAR_CONCENTRATION_UNITS,
    MASS_CONCENTRATION_UNITS
)


class TestVolumeConversion(unittest.TestCase):
    """Test volume unit conversions."""
    
    def test_ml_to_ul(self):
        """Test milliliter to microliter conversion."""
        result = convert_volume(1, 'mL', 'µL')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_ul_to_ml(self):
        """Test microliter to milliliter conversion."""
        result = convert_volume(1000, 'µL', 'mL')
        self.assertAlmostEqual(result, 1, places=6)
    
    def test_l_to_ml(self):
        """Test liter to milliliter conversion."""
        result = convert_volume(1, 'L', 'mL')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_nl_to_ul(self):
        """Test nanoliter to microliter conversion."""
        result = convert_volume(1000, 'nL', 'µL')
        self.assertAlmostEqual(result, 1, places=6)
    
    def test_same_unit_conversion(self):
        """Test conversion to same unit."""
        result = convert_volume(42, 'mL', 'mL')
        self.assertAlmostEqual(result, 42, places=6)
    
    def test_alternative_ul_spelling(self):
        """Test alternative microliter spelling (uL)."""
        result = convert_volume(1, 'mL', 'uL')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_invalid_unit(self):
        """Test that invalid units raise ValueError."""
        with self.assertRaises(ValueError):
            convert_volume(1, 'mL', 'invalid')
        with self.assertRaises(ValueError):
            convert_volume(1, 'invalid', 'mL')


class TestConcentrationConversion(unittest.TestCase):
    """Test concentration unit conversions."""
    
    def test_m_to_mm(self):
        """Test Molar to millimolar conversion."""
        result = convert_concentration(1, 'M', 'mM')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_mm_to_um(self):
        """Test millimolar to micromolar conversion."""
        result = convert_concentration(1, 'mM', 'µM')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_um_to_nm(self):
        """Test micromolar to nanomolar conversion."""
        result = convert_concentration(1, 'µM', 'nM')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_nm_to_m(self):
        """Test nanomolar to Molar conversion."""
        result = convert_concentration(1e9, 'nM', 'M')
        self.assertAlmostEqual(result, 1, places=6)
    
    def test_same_unit_conversion(self):
        """Test conversion to same unit."""
        result = convert_concentration(5.5, 'mM', 'mM')
        self.assertAlmostEqual(result, 5.5, places=6)
    
    def test_alternative_um_spelling(self):
        """Test alternative micromolar spelling (uM)."""
        result = convert_concentration(1, 'mM', 'uM')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_invalid_unit(self):
        """Test that invalid units raise ValueError."""
        with self.assertRaises(ValueError):
            convert_concentration(1, 'M', 'invalid')


class TestMassConcentrationConversion(unittest.TestCase):
    """Test mass-based concentration unit conversions."""
    
    def test_mg_ml_to_ug_ml(self):
        """Test mg/mL to µg/mL conversion."""
        result = convert_concentration(1, 'mg/mL', 'µg/mL')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_g_ml_to_mg_ml(self):
        """Test g/mL to mg/mL conversion."""
        result = convert_concentration(1, 'g/mL', 'mg/mL')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_ug_ml_to_ng_ml(self):
        """Test µg/mL to ng/mL conversion."""
        result = convert_concentration(1, 'µg/mL', 'ng/mL')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_g_l_to_mg_ml(self):
        """Test g/L to mg/mL conversion."""
        result = convert_concentration(1, 'g/L', 'mg/mL')
        self.assertAlmostEqual(result, 1, places=6)
    
    def test_mg_l_to_ug_ml(self):
        """Test mg/L to µg/mL conversion."""
        result = convert_concentration(1, 'mg/L', 'µg/mL')
        self.assertAlmostEqual(result, 1, places=6)
    
    def test_alternative_ug_spelling(self):
        """Test alternative microgram spelling (ug/mL)."""
        result = convert_concentration(1, 'mg/mL', 'ug/mL')
        self.assertAlmostEqual(result, 1000, places=6)
    
    def test_molar_to_mass_raises_error(self):
        """Test that converting between molar and mass units raises error."""
        with self.assertRaises(ValueError):
            convert_concentration(1, 'M', 'mg/mL')
        with self.assertRaises(ValueError):
            convert_concentration(1, 'mg/mL', 'µM')


class TestBasicDilution(unittest.TestCase):
    """Test basic dilution calculations (C1V1 = C2V2)."""
    
    def test_simple_dilution(self):
        """Test a simple 10x dilution."""
        v1, added = calculate_dilution(10, 1, 100)
        self.assertAlmostEqual(v1, 10, places=6)
        self.assertAlmostEqual(added, 90, places=6)
    
    def test_2x_dilution(self):
        """Test a 2x dilution."""
        v1, added = calculate_dilution(2, 1, 100)
        self.assertAlmostEqual(v1, 50, places=6)
        self.assertAlmostEqual(added, 50, places=6)
    
    def test_1000x_dilution(self):
        """Test a 1000x dilution."""
        v1, added = calculate_dilution(1000, 1, 1000)
        self.assertAlmostEqual(v1, 1, places=6)
        self.assertAlmostEqual(added, 999, places=6)
    
    def test_no_dilution(self):
        """Test when C1 equals C2 (no dilution needed)."""
        v1, added = calculate_dilution(5, 5, 100)
        self.assertAlmostEqual(v1, 100, places=6)
        self.assertAlmostEqual(added, 0, places=6)
    
    def test_fractional_concentrations(self):
        """Test with fractional concentration values."""
        v1, added = calculate_dilution(0.5, 0.1, 50)
        self.assertAlmostEqual(v1, 10, places=6)
        self.assertAlmostEqual(added, 40, places=6)
    
    def test_zero_c1_raises_error(self):
        """Test that C1=0 raises ValueError."""
        with self.assertRaises(ValueError) as context:
            calculate_dilution(0, 1, 100)
        self.assertIn("cannot be zero", str(context.exception))


class TestDilutionWithUnits(unittest.TestCase):
    """Test dilution calculations with different units."""
    
    def test_dilution_with_mm_units(self):
        """Test dilution with millimolar concentrations."""
        # 100 mM -> 10 mM in 500 µL
        v1, added = calculate_dilution(100, 10, 500, 'mM', 'mM', 'µL', 'µL')
        self.assertAlmostEqual(v1, 50, places=6)
        self.assertAlmostEqual(added, 450, places=6)
    
    def test_dilution_mixed_concentration_units(self):
        """Test dilution with different concentration units."""
        # 1 M -> 10 mM in 100 mL
        v1, added = calculate_dilution(1, 10, 100, 'M', 'mM', 'mL', 'mL')
        self.assertAlmostEqual(v1, 1, places=6)
        self.assertAlmostEqual(added, 99, places=6)
    
    def test_dilution_mixed_volume_units(self):
        """Test dilution with different volume units."""
        # Output in µL when input is mL
        v1, added = calculate_dilution(10, 1, 1, 'M', 'M', 'mL', 'µL')
        self.assertAlmostEqual(v1, 100, places=6)
        self.assertAlmostEqual(added, 900, places=6)
    
    def test_dilution_um_to_nm(self):
        """Test dilution from µM to nM."""
        # 100 µM -> 500 nM in 1 mL
        v1, added = calculate_dilution(100, 500, 1, 'µM', 'nM', 'mL', 'µL')
        self.assertAlmostEqual(v1, 5, places=6)
        self.assertAlmostEqual(added, 995, places=6)
    
    def test_dilution_output_in_different_unit(self):
        """Test that output unit conversion works correctly."""
        # Calculate in mL but output in L
        v1, added = calculate_dilution(10, 1, 100, 'M', 'M', 'mL', 'L')
        self.assertAlmostEqual(v1, 0.01, places=6)
        self.assertAlmostEqual(added, 0.09, places=6)
    
    def test_dilution_with_mass_units(self):
        """Test dilution using mass concentration units."""
        # 50 mg/mL -> 5 mg/mL in 100 mL
        v1, added = calculate_dilution(50, 5, 100, 'mg/mL', 'mg/mL', 'mL', 'mL')
        self.assertAlmostEqual(v1, 10, places=6)
        self.assertAlmostEqual(added, 90, places=6)
    
    def test_dilution_mass_mixed_units(self):
        """Test dilution with different mass units."""
        # 1 g/mL -> 100 mg/mL in 500 µL
        v1, added = calculate_dilution(1, 100, 500, 'g/mL', 'mg/mL', 'µL', 'µL')
        self.assertAlmostEqual(v1, 50, places=6)
        self.assertAlmostEqual(added, 450, places=6)
    
    def test_dilution_ug_ml_to_ng_ml(self):
        """Test dilution from µg/mL to ng/mL."""
        # 500 µg/mL -> 50 ng/mL in 10 mL (10000-fold dilution)
        v1, added = calculate_dilution(500, 50, 10, 'µg/mL', 'ng/mL', 'mL', 'mL')
        self.assertAlmostEqual(v1, 0.001, places=6)
        self.assertAlmostEqual(added, 9.999, places=3)
    
    def test_mixing_molar_and_mass_raises_error(self):
        """Test that mixing molar and mass units raises ValueError."""
        with self.assertRaises(ValueError):
            calculate_dilution(10, 5, 100, 'M', 'mg/mL', 'mL', 'mL')
        with self.assertRaises(ValueError):
            calculate_dilution(50, 10, 100, 'mg/mL', 'mM', 'mL', 'mL')


class TestValidateDilution(unittest.TestCase):
    """Test the dilution validation function."""
    
    def test_valid_dilution(self):
        """Test validation of a correct dilution."""
        # 10M * 10mL = 1M * 100mL
        self.assertTrue(validate_dilution(10, 1, 10, 100))
    
    def test_invalid_dilution(self):
        """Test validation of an incorrect dilution."""
        # 10M * 20mL ≠ 1M * 100mL
        self.assertFalse(validate_dilution(10, 1, 20, 100))
    
    def test_validation_with_tolerance(self):
        """Test validation with custom tolerance."""
        # Slight error within tolerance
        self.assertTrue(validate_dilution(10, 1, 10.05, 100, tolerance=0.01))
        # Error exceeds tolerance
        self.assertFalse(validate_dilution(10, 1, 12, 100, tolerance=0.01))
    
    def test_validation_zero_values(self):
        """Test validation with zero values."""
        self.assertTrue(validate_dilution(0, 0, 10, 100))
        self.assertTrue(validate_dilution(1, 0, 0, 100))


class TestRealWorldScenarios(unittest.TestCase):
    """Test real-world lab scenarios."""
    
    def test_protein_dilution(self):
        """Test typical protein stock dilution."""
        # Dilute 10 mg/mL protein to 1 mg/mL in 5 mL
        # (treating as 10 M to 1 M for the calculation)
        v1, added = calculate_dilution(10, 1, 5, 'M', 'M', 'mL', 'mL')
        self.assertAlmostEqual(v1, 0.5, places=6)
        self.assertAlmostEqual(added, 4.5, places=6)
    
    def test_pcr_master_mix(self):
        """Test PCR master mix dilution."""
        # 10X master mix to 1X in 50 µL
        v1, added = calculate_dilution(10, 1, 50, 'M', 'M', 'µL', 'µL')
        self.assertAlmostEqual(v1, 5, places=6)
        self.assertAlmostEqual(added, 45, places=6)
    
    def test_antibody_dilution(self):
        """Test antibody dilution (1:1000)."""
        # 1000X stock to 1X in 10 mL
        v1, added = calculate_dilution(1000, 1, 10, 'M', 'M', 'mL', 'µL')
        self.assertAlmostEqual(v1, 10, places=6)
        self.assertAlmostEqual(added, 9990, places=6)
    
    def test_serial_dilution_step(self):
        """Test one step of a serial dilution."""
        # 100 µM to 10 µM in 1 mL (10-fold dilution)
        v1, added = calculate_dilution(100, 10, 1, 'µM', 'µM', 'mL', 'µL')
        self.assertAlmostEqual(v1, 100, places=6)
        self.assertAlmostEqual(added, 900, places=6)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and potential error conditions."""
    
    def test_very_small_volumes(self):
        """Test calculations with very small volumes."""
        v1, added = calculate_dilution(1000, 1, 1, 'M', 'M', 'µL', 'nL')
        self.assertAlmostEqual(v1, 1, places=6)
        self.assertAlmostEqual(added, 999, places=6)
    
    def test_very_large_volumes(self):
        """Test calculations with large volumes."""
        v1, added = calculate_dilution(100, 1, 1, 'M', 'M', 'L', 'mL')
        self.assertAlmostEqual(v1, 10, places=6)
        self.assertAlmostEqual(added, 990, places=6)
    
    def test_high_concentration_ratios(self):
        """Test with very high concentration ratios."""
        v1, added = calculate_dilution(10000, 1, 1000, 'M', 'M', 'mL', 'µL')
        self.assertAlmostEqual(v1, 100, places=6)
        self.assertAlmostEqual(added, 999900, places=6)
    
    def test_c2_greater_than_c1_error(self):
        """Test that C2 > C1 gives negative added volume (concentration, not dilution)."""
        v1, added = calculate_dilution(1, 10, 100)
        # This would require concentrating, not diluting
        self.assertGreater(v1, 100)
        self.assertLess(added, 0)


def run_tests():
    """Run all tests and display results."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add all test classes
    suite.addTests(loader.loadTestsFromTestCase(TestVolumeConversion))
    suite.addTests(loader.loadTestsFromTestCase(TestConcentrationConversion))
    suite.addTests(loader.loadTestsFromTestCase(TestMassConcentrationConversion))
    suite.addTests(loader.loadTestsFromTestCase(TestBasicDilution))
    suite.addTests(loader.loadTestsFromTestCase(TestDilutionWithUnits))
    suite.addTests(loader.loadTestsFromTestCase(TestValidateDilution))
    suite.addTests(loader.loadTestsFromTestCase(TestRealWorldScenarios))
    suite.addTests(loader.loadTestsFromTestCase(TestEdgeCases))
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Return exit code
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    sys.exit(run_tests())
