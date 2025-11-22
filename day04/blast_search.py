"""
BLAST Search Module
This module provides functionality to search DNA sequences using NCBI BLAST API
"""

from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO
import time
import ssl
import urllib.request
import certifi
from datetime import datetime
import os


class BLASTSearcher:
    """Class to handle BLAST sequence searches via NCBI"""
    
    def __init__(self):
        self.database = "nt"  # nucleotide database
        self.program = "blastn"  # nucleotide BLAST
        self._setup_ssl()
    
    def _setup_ssl(self):
        """Setup SSL context to handle certificate verification"""
        try:
            # Try to use certifi's certificate bundle
            ssl_context = ssl.create_default_context(cafile=certifi.where())
            urllib.request.install_opener(
                urllib.request.build_opener(
                    urllib.request.HTTPSHandler(context=ssl_context)
                )
            )
        except Exception:
            # Fallback: create unverified context (not recommended for production)
            ssl._create_default_https_context = ssl._create_unverified_context
        
    def search_sequence(self, sequence, organism=None, max_results=10):
        """
        Perform BLAST search for a DNA sequence
        
        Args:
            sequence (str): DNA sequence to search
            organism (str, optional): Organism name to limit search
            max_results (int): Maximum number of results to return
            
        Returns:
            list: List of dictionaries containing alignment results
        """
        try:
            # Prepare entrez query for organism filtering
            entrez_query = None
            if organism and organism.strip():
                # Format: "Homo sapiens"[Organism]
                entrez_query = f'"{organism}"[Organism]'
            
            print(f"Submitting BLAST query...")
            print(f"Sequence length: {len(sequence)} bp")
            if entrez_query:
                print(f"Organism filter: {organism}")
            
            # Submit BLAST query
            result_handle = NCBIWWW.qblast(
                program=self.program,
                database=self.database,
                sequence=sequence,
                entrez_query=entrez_query,
                hitlist_size=max_results
            )
            
            # Parse results
            blast_records = NCBIXML.read(result_handle)
            result_handle.close()
            
            # Extract alignment information
            results = []
            for alignment in blast_records.alignments[:max_results]:
                for hsp in alignment.hsps:
                    result = {
                        'title': alignment.title,
                        'accession': alignment.accession,
                        'length': alignment.length,
                        'e_value': hsp.expect,
                        'score': hsp.score,
                        'identities': hsp.identities,
                        'positives': hsp.positives,
                        'gaps': hsp.gaps,
                        'alignment_length': hsp.align_length,
                        'query': hsp.query,
                        'match': hsp.match,
                        'subject': hsp.sbjct,
                        'identity_percent': (hsp.identities / hsp.align_length) * 100
                    }
                    results.append(result)
                    break  # Only take the first HSP per alignment
            
            return results
            
        except Exception as e:
            raise Exception(f"BLAST search failed: {str(e)}")
    
    def validate_sequence(self, sequence):
        """
        Validate DNA sequence
        
        Args:
            sequence (str): DNA sequence to validate
            
        Returns:
            tuple: (is_valid, error_message)
        """
        if not sequence or not sequence.strip():
            return False, "Sequence cannot be empty"
        
        # Remove whitespace
        clean_seq = ''.join(sequence.split()).upper()
        
        # Check for valid DNA characters
        valid_chars = set('ATCGN')
        invalid_chars = set(clean_seq) - valid_chars
        
        if invalid_chars:
            return False, f"Invalid DNA characters found: {', '.join(sorted(invalid_chars))}"
        
        if len(clean_seq) < 10:
            return False, "Sequence too short (minimum 10 nucleotides)"
        
        return True, clean_seq
    
    def format_result_summary(self, result):
        """
        Format a single result for display
        
        Args:
            result (dict): Result dictionary
            
        Returns:
            str: Formatted result string
        """
        summary = f"Title: {result['title'][:100]}...\n"
        summary += f"Accession: {result['accession']}\n"
        summary += f"Identity: {result['identity_percent']:.2f}%\n"
        summary += f"E-value: {result['e_value']:.2e}\n"
        summary += f"Score: {result['score']}\n"
        summary += f"Alignment Length: {result['alignment_length']}\n"
        summary += "-" * 50 + "\n"
        return summary
    
    def save_top_sequences_as_fasta(self, results, output_dir=".", top_n=3):
        """
        Save the top N sequences as FASTA files based on score
        
        Args:
            results (list): List of result dictionaries
            output_dir (str): Directory to save FASTA files
            top_n (int): Number of top sequences to save (default: 3)
            
        Returns:
            list: List of file paths that were created
        """
        if not results:
            return []
        
        # Sort results by score (descending)
        sorted_results = sorted(results, key=lambda x: x['score'], reverse=True)
        
        # Take top N results
        top_results = sorted_results[:top_n]
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate timestamp for unique filenames
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        saved_files = []
        
        for i, result in enumerate(top_results, 1):
            # Create SeqRecord object
            seq_record = SeqRecord(
                Seq(result['subject'].replace('-', '')),  # Remove gaps from alignment
                id=result['accession'],
                name=result['accession'],
                description=f"Rank {i} | Score: {result['score']} | Identity: {result['identity_percent']:.2f}% | {result['title']}"
            )
            
            # Generate filename
            safe_accession = result['accession'].replace('|', '_')
            filename = f"top{i}_{safe_accession}_{timestamp}.fasta"
            filepath = os.path.join(output_dir, filename)
            
            # Write FASTA file
            with open(filepath, 'w') as fasta_file:
                SeqIO.write(seq_record, fasta_file, 'fasta')
            
            saved_files.append(filepath)
            print(f"Saved: {filepath}")
        
        return saved_files


if __name__ == "__main__":
    # Test example
    searcher = BLASTSearcher()
    
    # Example DNA sequence (human insulin gene fragment)
    test_sequence = "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAGACGCAGCCCGCAGGCAGCCCCACACCCGCCGCCTCCTGCACCGAGAGAGATGGAATAAAGCCCTTGAACCAGC"
    
    print("Testing BLAST search...")
    is_valid, result = searcher.validate_sequence(test_sequence)
    
    if is_valid:
        print(f"Searching for sequence...")
        results = searcher.search_sequence(result, organism="Homo sapiens", max_results=5)
        print(f"\nFound {len(results)} results:\n")
        
        for i, res in enumerate(results[:3], 1):
            print(f"Result {i}:")
            print(searcher.format_result_summary(res))
    else:
        print(f"Invalid sequence: {result}")
