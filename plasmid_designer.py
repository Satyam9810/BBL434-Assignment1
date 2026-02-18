#!/usr/bin/env python3
"""
Plasmid Design Tool
===================
1. Input DNA sequence (finds ORI)
2. Design specifications (MCS and markers)
3. Marker database
"""

import re
import sys
from typing import Dict, List, Tuple, Optional
from collections import defaultdict


class PlasmidDesigner:
    """Main class for designing plasmids"""
    
    def __init__(self, markers_file: str):
        """
        Initialize the plasmid designer
        
        Args:
            markers_file: Path to the markers database file
        """
        self.restriction_sites = {}
        self.markers = self.load_markers(markers_file)
        
    def load_markers(self, markers_file: str) -> Dict[str, Dict[str, str]]:
        """
        Load markers from the database file
        
        Args:
            markers_file: Path to markers.tab file
            
        Returns:
            Dictionary of markers with their sequences and descriptions
        """
        markers = {}
        
        try:
            with open(markers_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 3:
                        name = parts[0]
                        marker_type = parts[1]
                        sequence = parts[2]
                        description = parts[3] if len(parts) > 3 else ""
                        
                        markers[name] = {
                            'type': marker_type,
                            'sequence': sequence,
                            'description': description
                        }
                        
                        # Store restriction enzyme recognition sites separately
                        if marker_type == 'RE':
                            self.restriction_sites[name] = sequence
                            
            print(f"Loaded {len(markers)} markers from database")
            print(f"  - {len(self.restriction_sites)} restriction enzymes")
            print(f"  - {len([m for m in markers.values() if m['type'] == 'AmpR'])} antibiotic markers")
            
        except FileNotFoundError:
            print(f"Warning: Markers file '{markers_file}' not found. Using default markers.")
            self._load_default_markers()
            
        return markers
    
    def _load_default_markers(self):
        # Basic restriction sites
        default_re = {
            'EcoRI': 'GAATTC',
            'BamHI': 'GGATCC',
            'HindIII': 'AAGCTT',
            'PstI': 'CTGCAG',
            'SalI': 'GTCGAC',
            'XbaI': 'TCTAGA',
            'KpnI': 'GGTACC',
            'SacI': 'GAGCTC',
            'SmaI': 'CCCGGG',
            'SphI': 'GCATGC'
        }
        
        for name, seq in default_re.items():
            self.markers[name] = {
                'type': 'RE',
                'sequence': seq,
                'description': f'{name} restriction site'
            }
            self.restriction_sites[name] = seq
    
    def read_fasta(self, fasta_file: str) -> Tuple[str, str]:
        """
        Read a FASTA file and return header and sequence
        
        Args:
            fasta_file: Path to FASTA file
            
        Returns:
            Tuple of (header, sequence)
        """
        header = ""
        sequence = ""
        
        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    header = line[1:]
                else:
                    sequence += line.upper()
                    
        return header, sequence
    
    def find_ori_by_pattern(self, sequence: str) -> Optional[str]:
        """
        Find origin of replication by looking for known ORI patterns
        
        Args:
            sequence: DNA sequence to search
            
        Returns:
            ORI sequence if found, None otherwise
        """
        # Look for pMB1/ColE1 origin signature
        # These origins have characteristic features like RNA II/RNA I regions
        
        # Pattern 1: Look for regions with high AT content (typical of ORI)
        window_size = 200
        max_at = 0
        best_ori_pos = 0
        
        for i in range(len(sequence) - window_size):
            window = sequence[i:i+window_size]
            at_content = (window.count('A') + window.count('T')) / len(window)
            
            if at_content > max_at:
                max_at = at_content
                best_ori_pos = i
        
        # Extract ~500bp region around the high AT content region
        ori_start = max(0, best_ori_pos - 150)
        ori_end = min(len(sequence), best_ori_pos + window_size + 150)
        ori_candidate = sequence[ori_start:ori_end]
        
        print(f"Found ORI candidate at position {ori_start}-{ori_end}")
        print(f"  AT content: {max_at*100:.1f}%")
        
        return ori_candidate
    
    def find_ori_by_gc_skew(self, sequence: str, window: int = 100) -> Optional[str]:
        """
        Find ORI using GC skew analysis
        GC skew = (G-C)/(G+C)
        ORI regions typically show characteristic GC skew patterns
        
        Args:
            sequence: DNA sequence
            window: Window size for calculation
            
        Returns:
            ORI sequence if found
        """
        gc_skew = []
        
        for i in range(0, len(sequence) - window, window):
            subseq = sequence[i:i+window]
            g_count = subseq.count('G')
            c_count = subseq.count('C')
            
            if (g_count + c_count) > 0:
                skew = (g_count - c_count) / (g_count + c_count)
            else:
                skew = 0
                
            gc_skew.append((i, skew))
        
        # Find minimum skew (characteristic of ORI)
        if gc_skew:
            min_pos = min(gc_skew, key=lambda x: x[1])[0]
            
            # Extract region around minimum
            ori_start = max(0, min_pos - 250)
            ori_end = min(len(sequence), min_pos + 250)
            
            print(f"GC skew minimum at position {min_pos}")
            return sequence[ori_start:ori_end]
        
        return None
    
    def find_ori_in_sequence(self, sequence: str) -> str:
        """
        Find origin of replication using multiple methods
        
        Args:
            sequence: DNA sequence to analyze
            
        Returns:
            ORI sequence
        """
        print("\n=== Searching for Origin of Replication (ORI) ===")
        
        # Method 1: Check if we have known ORI in markers
        if 'ori_pMB1' in self.markers:
            ori_seq = self.markers['ori_pMB1']['sequence']
            if ori_seq in sequence:
                print("Found known ori_pMB1 in sequence!")
                # Extract with some flanking region
                pos = sequence.find(ori_seq)
                start = max(0, pos - 50)
                end = min(len(sequence), pos + len(ori_seq) + 50)
                return sequence[start:end]
        
        # Method 2: Try pattern-based search
        ori_pattern = self.find_ori_by_pattern(sequence)
        
        # Method 3: Try GC skew method
        ori_gc_skew = self.find_ori_by_gc_skew(sequence)
        
        # Use the pattern-based result if available
        if ori_pattern:
            return ori_pattern
        elif ori_gc_skew:
            return ori_gc_skew
        else:
            # Fallback: use first 500bp as ORI region
            print("Warning: Could not definitively identify ORI. Using first 500bp.")
            return sequence[:500]
    
    def parse_design_file(self, design_file: str) -> Dict[str, List[Tuple[str, str]]]:
        """
        Parse the design file to extract plasmid components
        
        Args:
            design_file: Path to design file
            
        Returns:
            Dictionary with 'mcs' and 'markers' lists
        """
        design = {
            'mcs': [],  # Multiple cloning sites
            'markers': [],  # Antibiotic resistance and other markers
            'origins': []  # Replication origins
        }
        
        with open(design_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                if ',' in line:
                    parts = [p.strip() for p in line.split(',')]
                    if len(parts) >= 2:
                        name = parts[0]
                        marker_name = parts[1]
                        
                        # Determine type based on name
                        if 'site' in name.lower() or marker_name in self.restriction_sites:
                            design['mcs'].append((name, marker_name))
                        elif 'ori' in name.lower() or 'origin' in marker_name.lower() or 'replication' in marker_name.lower():
                            design['origins'].append((name, marker_name))
                        else:
                            design['markers'].append((name, marker_name))
        
        print(f"\nParsed design file:")
        print(f"  - {len(design['mcs'])} restriction sites")
        print(f"  - {len(design['markers'])} markers")
        print(f"  - {len(design['origins'])} origins")
        
        return design
    
    def build_mcs(self, mcs_sites: List[Tuple[str, str]]) -> str:
        """
        Build the multiple cloning site region
        
        Args:
            mcs_sites: List of (name, enzyme) tuples
            
        Returns:
            MCS sequence
        """
        mcs_sequence = ""
        
        # Add spacer before MCS
        mcs_sequence += "GGGCCC"  # SmaI site as spacer
        
        for name, enzyme in mcs_sites:
            if enzyme in self.markers:
                site_seq = self.markers[enzyme]['sequence']
                mcs_sequence += site_seq
                print(f"  Added {enzyme} site: {site_seq}")
            elif enzyme in self.restriction_sites:
                site_seq = self.restriction_sites[enzyme]
                mcs_sequence += site_seq
                print(f"  Added {enzyme} site: {site_seq}")
            else:
                print(f"  Warning: {enzyme} not found in markers database")
        
        # Add spacer after MCS
        mcs_sequence += "GGGCCC"
        
        return mcs_sequence
    
    def construct_plasmid(self, ori_seq: str, design: Dict, 
                         remove_sites: Optional[List[str]] = None) -> str:
        """
        Construct the final plasmid sequence
        
        Args:
            ori_seq: Origin of replication sequence
            design: Design dictionary with components
            remove_sites: List of restriction sites to remove
            
        Returns:
            Complete plasmid sequence
        """
        print("\n=== Constructing Plasmid ===")
        
        plasmid_parts = []
        
        # 1. Add origin of replication
        print("Adding origin of replication...")
        plasmid_parts.append(("ORI", ori_seq))
        
        # 2. Add antibiotic resistance markers
        print("Adding antibiotic resistance markers...")
        for name, marker_name in design['markers']:
            # Try direct match first
            marker_key = marker_name
            if marker_key not in self.markers:
                # Try alternative keys
                for key in self.markers:
                    if marker_name.lower() in key.lower() or key.lower() in marker_name.lower():
                        marker_key = key
                        break
            
            if marker_key in self.markers:
                marker_seq = self.markers[marker_key]['sequence']
                plasmid_parts.append((name, marker_seq))
                print(f"  Added {marker_name}: {len(marker_seq)} bp")
            else:
                print(f"  Warning: Marker '{marker_name}' not found in database")
        
        # 3. Build and add MCS
        if design['mcs']:
            print("Building multiple cloning site...")
            mcs_seq = self.build_mcs(design['mcs'])
            plasmid_parts.append(("MCS", mcs_seq))
        
        # 4. Combine all parts
        plasmid_sequence = ""
        for part_name, part_seq in plasmid_parts:
            plasmid_sequence += part_seq
        
        # 5. Remove specified restriction sites from non-MCS regions
        if remove_sites:
            print(f"\nRemoving restriction sites: {', '.join(remove_sites)}")
            plasmid_sequence = self.remove_restriction_sites(
                plasmid_sequence, remove_sites, design['mcs']
            )
        
        print(f"\nFinal plasmid size: {len(plasmid_sequence)} bp")
        
        return plasmid_sequence
    
    def remove_restriction_sites(self, sequence: str, sites_to_remove: List[str],
                                mcs_sites: List[Tuple[str, str]]) -> str:
        """
        Remove specified restriction sites from sequence
        (but preserve them in MCS if they're part of the design)
        
        Args:
            sequence: DNA sequence
            sites_to_remove: List of restriction enzyme names to remove
            mcs_sites: List of MCS sites to preserve
            
        Returns:
            Modified sequence
        """
        # Get list of sites that should be preserved (in MCS)
        preserve_sites = [enzyme for name, enzyme in mcs_sites]
        
        modified_seq = sequence
        
        for site_name in sites_to_remove:
            # Skip if this site should be preserved in MCS
            if site_name in preserve_sites:
                continue
                
            if site_name in self.restriction_sites:
                site_seq = self.restriction_sites[site_name]
                
                # Count occurrences
                count = modified_seq.count(site_seq)
                
                if count > 0:
                    # Simple removal: introduce silent mutation
                    # For this implementation, we'll just remove the site
                    # In practice, you'd want to maintain reading frames
                    modified_seq = modified_seq.replace(site_seq, 'N' * len(site_seq))
                    print(f"  Removed {count} occurrence(s) of {site_name} ({site_seq})")
        
        return modified_seq
    
    def write_fasta(self, filename: str, header: str, sequence: str):
        """
        Write sequence to FASTA file
        
        Args:
            filename: Output filename
            header: FASTA header
            sequence: DNA sequence
        """
        with open(filename, 'w') as f:
            f.write(f">{header}\n")
            
            # Write sequence in 60-character lines
            for i in range(0, len(sequence), 60):
                f.write(sequence[i:i+60] + '\n')
        
        print(f"\nOutput written to: {filename}")
    
    def analyze_plasmid(self, sequence: str) -> Dict:
        """
        Analyze plasmid sequence for various features
        
        Args:
            sequence: Plasmid sequence
            
        Returns:
            Dictionary with analysis results
        """
        analysis = {
            'length': len(sequence),
            'gc_content': (sequence.count('G') + sequence.count('C')) / len(sequence) * 100,
            'restriction_sites': {}
        }
        
        # Find restriction sites
        for enzyme, site in self.restriction_sites.items():
            count = sequence.count(site)
            if count > 0:
                analysis['restriction_sites'][enzyme] = count
        
        return analysis
    
    def print_analysis(self, analysis: Dict):
        """Print plasmid analysis results"""
        print("\n=== Plasmid Analysis ===")
        print(f"Length: {analysis['length']} bp")
        print(f"GC content: {analysis['gc_content']:.1f}%")
        print(f"\nRestriction sites found:")
        
        if analysis['restriction_sites']:
            for enzyme, count in sorted(analysis['restriction_sites'].items()):
                print(f"  {enzyme}: {count} site(s)")
        else:
            print("  None found")


def main():
    """Main function"""
    print("="*60)
    print("Plasmid Design Tool")
    print("="*60)
    
    # Check command line arguments
    if len(sys.argv) < 4:
        print("\nUsage: python plasmid_designer.py <input.fa> <design.txt> <output.fa> [markers.tab]")
        print("\nExample:")
        print("  python plasmid_designer.py pUC19.fa Design_pUC19.txt Output.fa markers.tab")
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    design_file = sys.argv[2]
    output_fasta = sys.argv[3]
    markers_file = sys.argv[4] if len(sys.argv) > 4 else "markers.tab"
    
    # Initialize designer
    designer = PlasmidDesigner(markers_file)
    
    # Read input sequence
    print(f"\nReading input sequence from: {input_fasta}")
    header, sequence = designer.read_fasta(input_fasta)
    print(f"Input sequence: {header}")
    print(f"Length: {len(sequence)} bp")
    
    # Find ORI
    ori_sequence = designer.find_ori_in_sequence(sequence)
    
    # Parse design
    print(f"\nReading design from: {design_file}")
    design = designer.parse_design_file(design_file)
    
    # Determine which sites to remove (EcoRI if not in design)
    all_enzymes_in_design = [enzyme for name, enzyme in design['mcs']]
    sites_to_remove = ['EcoRI'] if 'EcoRI' not in all_enzymes_in_design else []
    
    # Construct plasmid
    plasmid_seq = designer.construct_plasmid(ori_sequence, design, sites_to_remove)
    
    # Analyze plasmid
    analysis = designer.analyze_plasmid(plasmid_seq)
    designer.print_analysis(analysis)
    
    # Write output
    output_header = f"Designed_Plasmid_from_{header}"
    designer.write_fasta(output_fasta, output_header, plasmid_seq)
    
    print("\n" + "="*60)
    print("Plasmid design completed successfully!")
    print("="*60)


if __name__ == "__main__":
    main()
