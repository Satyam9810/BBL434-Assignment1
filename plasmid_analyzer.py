#!/usr/bin/env python3
"""
Plasmid Analyzer
================

Analyzes a plasmid sequence and provides detailed statistics
and feature detection.
"""

import sys
from collections import Counter


class PlasmidAnalyzer:
    
    def __init__(self):
        self.restriction_enzymes = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'PstI': 'CTGCAG',
            'SalI': 'GTCGAC',
            'XbaI': 'TCTAGA',
            'KpnI': 'GGTACC',
            'SacI': 'GAGCTC',
            'SmaI': 'CCCGGG',
            'SphI': 'GCATGC',
            'NotI': 'GCGGCCGC',
            'XhoI': 'CTCGAG',
            'NcoI': 'CCATGG',
            'NdeI': 'CATATG',
            'BglII': 'AGATCT',
            'ApaI': 'GGGCCC',
            'SpeI': 'ACTAGT',
            'PvuII': 'CAGCTG'
        }
    
    def read_fasta(self, filename):
        header = ""
        sequence = ""
        
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    header = line[1:]
                else:
                    sequence += line.upper()
        
        return header, sequence
    
    def analyze_composition(self, sequence):
        counts = Counter(sequence)
        total = len(sequence)
        
        composition = {
            'A': counts.get('A', 0),
            'T': counts.get('T', 0),
            'G': counts.get('G', 0),
            'C': counts.get('C', 0),
            'N': counts.get('N', 0),
            'Other': total - sum([counts.get(n, 0) for n in 'ATGCN'])
        }
        
        # Calculate percentages
        percentages = {base: (count/total*100) for base, count in composition.items()}
        
        # GC content
        gc_content = (composition['G'] + composition['C']) / total * 100
        
        # AT content
        at_content = (composition['A'] + composition['T']) / total * 100
        
        return composition, percentages, gc_content, at_content
    
    def find_restriction_sites(self, sequence):
        sites = {}
        
        for enzyme, site_seq in self.restriction_enzymes.items():
            count = sequence.count(site_seq)
            if count > 0:
                # Find positions
                positions = []
                start = 0
                while True:
                    pos = sequence.find(site_seq, start)
                    if pos == -1:
                        break
                    positions.append(pos)
                    start = pos + 1
                
                sites[enzyme] = {
                    'count': count,
                    'sequence': site_seq,
                    'positions': positions
                }
        
        return sites
    
    def analyze_orfs(self, sequence):
        # Start codons
        start_codons = ['ATG']
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        orfs = []
        
        # Check all 3 reading frames
        for frame in range(3):
            pos = frame
            while pos < len(sequence) - 2:
                codon = sequence[pos:pos+3]
                
                if codon in start_codons:
                    # Find stop codon
                    orf_start = pos
                    orf_pos = pos + 3
                    
                    while orf_pos < len(sequence) - 2:
                        stop_codon = sequence[orf_pos:orf_pos+3]
                        if stop_codon in stop_codons:
                            orf_length = orf_pos - orf_start + 3
                            if orf_length >= 300:  # At least 100 amino acids
                                orfs.append({
                                    'start': orf_start,
                                    'end': orf_pos + 3,
                                    'length': orf_length,
                                    'frame': frame,
                                    'aa_length': orf_length // 3
                                })
                            break
                        orf_pos += 3
                
                pos += 3
        
        return orfs
    
    def print_analysis(self, header, sequence):
        """Print comprehensive analysis"""
        print("="*70)
        print("PLASMID ANALYSIS REPORT")
        print("="*70)
        print(f"\nPlasmid: {header}")
        print(f"Length: {len(sequence)} bp")
        
        # Composition analysis
        print("\n" + "-"*70)
        print("NUCLEOTIDE COMPOSITION")
        print("-"*70)
        
        comp, perc, gc, at = self.analyze_composition(sequence)
        
        print(f"\nBase counts:")
        for base, count in comp.items():
            if count > 0:
                print(f"  {base}: {count:6d} ({perc[base]:5.2f}%)")
        
        print(f"\nGC content: {gc:.2f}%")
        print(f"AT content: {at:.2f}%")
        
        # Restriction sites
        print("\n" + "-"*70)
        print("RESTRICTION ENZYME SITES")
        print("-"*70)
        
        sites = self.find_restriction_sites(sequence)
        
        if sites:
            print(f"\nFound {len(sites)} unique enzyme sites:")
            print(f"\n{'Enzyme':<12} {'Sites':<6} {'Sequence':<10} {'Positions'}")
            print("-"*70)
            
            for enzyme in sorted(sites.keys()):
                info = sites[enzyme]
                pos_str = ', '.join(str(p) for p in info['positions'][:5])
                if len(info['positions']) > 5:
                    pos_str += f" ... (+{len(info['positions'])-5} more)"
                
                print(f"{enzyme:<12} {info['count']:<6} {info['sequence']:<10} {pos_str}")
        else:
            print("\nNo restriction sites found for common enzymes")
        
        # ORF analysis
        print("\n" + "-"*70)
        print("OPEN READING FRAMES (ORFs)")
        print("-"*70)
        
        orfs = self.analyze_orfs(sequence)
        
        if orfs:
            print(f"\nFound {len(orfs)} ORFs (≥ 300 bp):")
            print(f"\n{'Start':<8} {'End':<8} {'Length':<10} {'Frame':<8} {'AA Length'}")
            print("-"*70)
            
            for orf in sorted(orfs, key=lambda x: x['length'], reverse=True)[:10]:
                print(f"{orf['start']:<8} {orf['end']:<8} {orf['length']:<10} "
                      f"{orf['frame']:<8} {orf['aa_length']}")
            
            if len(orfs) > 10:
                print(f"\n... and {len(orfs)-10} more ORFs")
        else:
            print("\nNo significant ORFs found (minimum 300 bp)")
        
        # Summary statistics
        print("\n" + "-"*70)
        print("SUMMARY STATISTICS")
        print("-"*70)
        
        print(f"\nTotal length:          {len(sequence)} bp")
        print(f"GC content:            {gc:.2f}%")
        print(f"Restriction sites:     {len(sites)}")
        print(f"ORFs (≥300bp):         {len(orfs)}")
        print(f"Ambiguous bases (N):   {comp['N']}")
        
        print("\n" + "="*70)


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 plasmid_analyzer.py <plasmid.fa>")
        print("\nExample:")
        print("  python3 plasmid_analyzer.py Output.fa")
        sys.exit(1)
    
    filename = sys.argv[1]
    
    analyzer = PlasmidAnalyzer()
    header, sequence = analyzer.read_fasta(filename)
    analyzer.print_analysis(header, sequence)


if __name__ == "__main__":
    main()
