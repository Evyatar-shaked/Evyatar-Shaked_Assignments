"""
Binary Gibson Assembly Wizard
Interactive tool for designing primers for binary Gibson assembly (vector + insert)
"""

import sys
from typing import Dict, List, Tuple
from primer_design import PrimerDesigner


class BinaryAssemblyWizard:
    """Wizard for binary Gibson assembly primer design"""
    
    def __init__(self):
        self.designer = PrimerDesigner()
        self.vector_sequence = ""
        self.insert_sequence = ""
        self.vector_cut_site = None
        self.primers = {}
        
    def print_header(self, text: str):
        """Print formatted header"""
        print("\n" + "=" * 70)
        print(text.center(70))
        print("=" * 70 + "\n")
    
    def print_section(self, text: str):
        """Print formatted section"""
        print("\n" + "-" * 70)
        print(text)
        print("-" * 70)
    
    def get_sequence_input(self, prompt: str) -> str:
        """Get and validate DNA sequence input"""
        while True:
            seq = input(prompt).strip().upper()
            # Remove whitespace and common formatting
            seq = seq.replace(" ", "").replace("\n", "").replace("\r", "")
            
            # Validate sequence
            valid_bases = set("ATGC")
            if all(base in valid_bases for base in seq):
                return seq
            else:
                print("❌ Invalid sequence. Please use only A, T, G, C characters.")
    
    def get_int_input(self, prompt: str, min_val: int = None, max_val: int = None) -> int:
        """Get and validate integer input"""
        while True:
            try:
                value = int(input(prompt).strip())
                if min_val is not None and value < min_val:
                    print(f"❌ Value must be at least {min_val}")
                    continue
                if max_val is not None and value > max_val:
                    print(f"❌ Value must be at most {max_val}")
                    continue
                return value
            except ValueError:
                print("❌ Please enter a valid integer.")
    
    def get_yes_no(self, prompt: str) -> bool:
        """Get yes/no input"""
        while True:
            response = input(prompt + " (y/n): ").strip().lower()
            if response in ['y', 'yes']:
                return True
            elif response in ['n', 'no']:
                return False
            else:
                print("❌ Please enter 'y' or 'n'")
    
    def display_validation_results(self, validation: Dict, primer_name: str):
        """Display validation results in a formatted way"""
        print(f"\n📊 Validation Results for {primer_name}:")
        print("-" * 50)
        
        metrics = validation['metrics']
        print(f"  Length: {metrics['length']} nt")
        print(f"  GC Content: {metrics['gc_content']}%")
        print(f"  Tm: {metrics['tm']}°C")
        print(f"  GC Clamp: {'✓ Yes' if metrics['has_gc_clamp'] else '✗ No'}")
        
        if 'annealing_region' in validation:
            anneal = validation['annealing_region']
            print(f"\n  Annealing Region:")
            print(f"    Length: {anneal['length']} nt")
            print(f"    Tm: {anneal['tm']}°C")
            print(f"    Sequence: {anneal['sequence']}")
        
        if 'homology_region' in validation:
            homology = validation['homology_region']
            print(f"\n  Homology Region:")
            print(f"    Length: {homology['length']} nt")
            print(f"    Tm: {homology['tm']}°C")
            print(f"    Sequence: {homology['sequence']}")
        
        # Display warnings
        if validation['warnings']:
            print(f"\n⚠️  Warnings:")
            for warning in validation['warnings']:
                print(f"    • {warning}")
        
        # Display errors
        if validation['errors']:
            print(f"\n❌ Errors:")
            for error in validation['errors']:
                print(f"    • {error}")
        
        if not validation['warnings'] and not validation['errors']:
            print("\n✅ All checks passed!")
    
    def run(self):
        """Run the interactive wizard"""
        self.print_header("🧬 Binary Gibson Assembly Primer Design Wizard")
        print("This wizard will help you design primers for binary Gibson assembly")
        print("following the Takara InFusion protocol (60°C incubation).\n")
        print("You will design primers for:")
        print("  1. Vector linearization (with homology for insert)")
        print("  2. Insert amplification (with homology for vector)")
        
        # Step 1: Get vector sequence
        self.print_section("Step 1: Vector Information")
        print("Enter your vector (plasmid) sequence.")
        print("This is the backbone that will receive the insert.\n")
        
        self.vector_sequence = self.get_sequence_input("Paste vector sequence: ")
        print(f"✓ Vector sequence loaded: {len(self.vector_sequence)} bp")
        
        # Step 2: Get linearization site
        self.print_section("Step 2: Vector Linearization Site")
        print("Where should the vector be linearized to accept the insert?")
        print("Enter the position (1-based) where the vector will be cut.\n")
        
        self.vector_cut_site = self.get_int_input(
            f"Cut site position (1-{len(self.vector_sequence)}): ",
            min_val=1,
            max_val=len(self.vector_sequence)
        ) - 1  # Convert to 0-based
        
        print(f"✓ Cut site set at position {self.vector_cut_site + 1}")
        
        # Step 3: Get insert sequence
        self.print_section("Step 3: Insert Information")
        print("Enter the insert (gene/fragment) sequence to be cloned into the vector.\n")
        
        self.insert_sequence = self.get_sequence_input("Paste insert sequence: ")
        print(f"✓ Insert sequence loaded: {len(self.insert_sequence)} bp")
        
        # Step 4: Design homology regions
        self.print_section("Step 4: Homology Region Design")
        print("Gibson assembly requires homology regions between fragments.")
        print(f"Recommended length: {self.designer.recommended_homology_length[0]}-{self.designer.recommended_homology_length[1]} nt")
        
        homology_length = self.get_int_input(
            f"\nHomology length ({self.designer.min_homology_length}-{self.designer.max_homology_length} nt) [default: 25]: ",
            min_val=self.designer.min_homology_length,
            max_val=self.designer.max_homology_length
        ) if self.get_yes_no("Customize homology length?") else 25
        
        # Design primers for vector linearization
        self.print_section("Step 5: Designing Vector Linearization Primers")
        print("These primers will linearize your vector and add homology for the insert.\n")
        
        # Forward primer: adds homology to insert start
        insert_start_homology = self.insert_sequence[:homology_length]
        
        # Reverse primer: adds homology to insert end
        insert_end_homology = self.insert_sequence[-homology_length:]
        
        print("Designing forward primer (upstream of cut site)...")
        vector_fwd_primer, vector_fwd_validation = self.designer.design_gibson_primer(
            self.vector_sequence,
            max(0, self.vector_cut_site - 50),  # Start 50 bp upstream
            insert_start_homology,
            direction='forward'
        )
        
        print("Designing reverse primer (downstream of cut site)...")
        vector_rev_primer, vector_rev_validation = self.designer.design_gibson_primer(
            self.vector_sequence,
            min(len(self.vector_sequence) - 1, self.vector_cut_site + 50),  # Start 50 bp downstream
            insert_end_homology,
            direction='reverse'
        )
        
        self.primers['vector_forward'] = vector_fwd_primer
        self.primers['vector_reverse'] = vector_rev_primer
        
        # Check Tm compatibility
        vector_fwd_tm = vector_fwd_validation['annealing_region']['tm']
        vector_rev_tm = vector_rev_validation['annealing_region']['tm']
        tm_diff = abs(vector_fwd_tm - vector_rev_tm)
        
        print(f"\n✓ Vector primers designed!")
        print(f"  Tm difference: {tm_diff:.1f}°C", end="")
        if tm_diff <= self.designer.max_tm_difference:
            print(" ✅ (within recommended 3°C)")
        else:
            print(f" ⚠️  (exceeds recommended {self.designer.max_tm_difference}°C)")
        
        # Design primers for insert amplification
        self.print_section("Step 6: Designing Insert Amplification Primers")
        print("These primers will amplify your insert and add homology for the vector.\n")
        
        # Get vector sequences at cut site for homology
        vector_upstream_homology = self.vector_sequence[
            max(0, self.vector_cut_site - homology_length):self.vector_cut_site
        ]
        vector_downstream_homology = self.vector_sequence[
            self.vector_cut_site:min(len(self.vector_sequence), self.vector_cut_site + homology_length)
        ]
        
        # If cut site is at edge, adjust homology
        if len(vector_upstream_homology) < homology_length:
            # Use from end of vector (circular)
            needed = homology_length - len(vector_upstream_homology)
            vector_upstream_homology = self.vector_sequence[-needed:] + vector_upstream_homology
        
        if len(vector_downstream_homology) < homology_length:
            # Use from start of vector (circular)
            needed = homology_length - len(vector_downstream_homology)
            vector_downstream_homology = vector_downstream_homology + self.vector_sequence[:needed]
        
        print("Designing forward primer (insert start)...")
        insert_fwd_primer, insert_fwd_validation = self.designer.design_gibson_primer(
            self.insert_sequence,
            0,  # Start of insert
            vector_upstream_homology,
            direction='forward'
        )
        
        print("Designing reverse primer (insert end)...")
        insert_rev_primer, insert_rev_validation = self.designer.design_gibson_primer(
            self.insert_sequence,
            len(self.insert_sequence) - 1,  # End of insert
            vector_downstream_homology,
            direction='reverse'
        )
        
        self.primers['insert_forward'] = insert_fwd_primer
        self.primers['insert_reverse'] = insert_rev_primer
        
        # Check Tm compatibility
        insert_fwd_tm = insert_fwd_validation['annealing_region']['tm']
        insert_rev_tm = insert_rev_validation['annealing_region']['tm']
        tm_diff = abs(insert_fwd_tm - insert_rev_tm)
        
        print(f"\n✓ Insert primers designed!")
        print(f"  Tm difference: {tm_diff:.1f}°C", end="")
        if tm_diff <= self.designer.max_tm_difference:
            print(" ✅ (within recommended 3°C)")
        else:
            print(f" ⚠️  (exceeds recommended {self.designer.max_tm_difference}°C)")
        
        # Display all results
        self.print_section("Step 7: Results Summary")
        
        print("\n🔬 VECTOR LINEARIZATION PRIMERS")
        print("=" * 70)
        print("\n1️⃣  Vector Forward Primer (VF):")
        print(f"   5'- {vector_fwd_primer} -3'")
        print(f"   Length: {len(vector_fwd_primer)} nt")
        self.display_validation_results(vector_fwd_validation, "Vector Forward")
        
        print("\n2️⃣  Vector Reverse Primer (VR):")
        print(f"   5'- {vector_rev_primer} -3'")
        print(f"   Length: {len(vector_rev_primer)} nt")
        self.display_validation_results(vector_rev_validation, "Vector Reverse")
        
        print("\n\n🧬 INSERT AMPLIFICATION PRIMERS")
        print("=" * 70)
        print("\n3️⃣  Insert Forward Primer (IF):")
        print(f"   5'- {insert_fwd_primer} -3'")
        print(f"   Length: {len(insert_fwd_primer)} nt")
        self.display_validation_results(insert_fwd_validation, "Insert Forward")
        
        print("\n4️⃣  Insert Reverse Primer (IR):")
        print(f"   5'- {insert_rev_primer} -3'")
        print(f"   Length: {len(insert_rev_primer)} nt")
        self.display_validation_results(insert_rev_validation, "Insert Reverse")
        
        # Assembly instructions
        self.print_section("Gibson Assembly Protocol")
        print("\n📋 PCR Steps:")
        print("  1. Amplify vector with VF + VR primers")
        print("  2. Amplify insert with IF + IR primers")
        print("  3. Verify products by gel electrophoresis")
        print("  4. Purify PCR products (gel extraction or column)")
        
        print("\n🧪 Gibson Assembly (Takara InFusion):")
        print("  1. Mix equimolar amounts of vector and insert PCR products")
        print("     - Recommended ratio: 1:1 to 1:2 (vector:insert)")
        print("     - Total DNA: 50-200 ng")
        print("  2. Add InFusion enzyme mix")
        print("  3. Incubate at 50°C for 15 minutes")
        print("  4. Transform into competent cells")
        
        print("\n💡 Tips:")
        print("  • Use high-fidelity polymerase for PCR")
        print("  • Gel-purify PCR products to remove template")
        print("  • Verify insert orientation by colony PCR or sequencing")
        print("  • DpnI digest can remove template plasmid")
        
        # Save to file option
        if self.get_yes_no("\nSave primers to file?"):
            self.save_primers_to_file(
                vector_fwd_primer, vector_rev_primer,
                insert_fwd_primer, insert_rev_primer,
                vector_fwd_validation, vector_rev_validation,
                insert_fwd_validation, insert_rev_validation
            )
        
        self.print_header("✅ Primer Design Complete!")
        print("Thank you for using the Binary Gibson Assembly Wizard!\n")
    
    def save_primers_to_file(self, vf, vr, inf, inr, vf_val, vr_val, inf_val, inr_val):
        """Save primers and validation results to a text file"""
        filename = input("Enter filename (e.g., primers.txt): ").strip()
        if not filename:
            filename = "gibson_primers.txt"
        
        try:
            with open(filename, 'w') as f:
                f.write("=" * 70 + "\n")
                f.write("Binary Gibson Assembly Primer Design Results\n")
                f.write("=" * 70 + "\n\n")
                
                f.write("VECTOR LINEARIZATION PRIMERS\n")
                f.write("-" * 70 + "\n")
                f.write(f"Vector Forward (VF):\n")
                f.write(f"  5'- {vf} -3'\n")
                f.write(f"  Length: {len(vf)} nt\n")
                f.write(f"  Tm (annealing): {vf_val['annealing_region']['tm']}°C\n")
                f.write(f"  GC%: {vf_val['metrics']['gc_content']}%\n\n")
                
                f.write(f"Vector Reverse (VR):\n")
                f.write(f"  5'- {vr} -3'\n")
                f.write(f"  Length: {len(vr)} nt\n")
                f.write(f"  Tm (annealing): {vr_val['annealing_region']['tm']}°C\n")
                f.write(f"  GC%: {vr_val['metrics']['gc_content']}%\n\n")
                
                f.write("\nINSERT AMPLIFICATION PRIMERS\n")
                f.write("-" * 70 + "\n")
                f.write(f"Insert Forward (IF):\n")
                f.write(f"  5'- {inf} -3'\n")
                f.write(f"  Length: {len(inf)} nt\n")
                f.write(f"  Tm (annealing): {inf_val['annealing_region']['tm']}°C\n")
                f.write(f"  GC%: {inf_val['metrics']['gc_content']}%\n\n")
                
                f.write(f"Insert Reverse (IR):\n")
                f.write(f"  5'- {inr} -3'\n")
                f.write(f"  Length: {len(inr)} nt\n")
                f.write(f"  Tm (annealing): {inr_val['annealing_region']['tm']}°C\n")
                f.write(f"  GC%: {inr_val['metrics']['gc_content']}%\n\n")
                
                f.write("\nPCR CONDITIONS\n")
                f.write("-" * 70 + "\n")
                f.write("Recommended annealing temperature:\n")
                avg_tm = (vf_val['annealing_region']['tm'] + vr_val['annealing_region']['tm']) / 2
                f.write(f"  Vector PCR: {avg_tm - 5:.1f}°C (Tm - 5°C)\n")
                avg_tm = (inf_val['annealing_region']['tm'] + inr_val['annealing_region']['tm']) / 2
                f.write(f"  Insert PCR: {avg_tm - 5:.1f}°C (Tm - 5°C)\n")
            
            print(f"✓ Primers saved to {filename}")
        except Exception as e:
            print(f"❌ Error saving file: {e}")


def main():
    """Main entry point"""
    wizard = BinaryAssemblyWizard()
    try:
        wizard.run()
    except KeyboardInterrupt:
        print("\n\n⚠️  Wizard interrupted by user.")
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ An error occurred: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
