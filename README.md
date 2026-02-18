# Plasmid Designer Tool

A comprehensive bioinformatics tool for designing plasmids based on genomic sequences and user specifications. This tool automatically detects origins of replication (ORI) and constructs functional plasmids with specified restriction sites and selectable markers.

## Features

- **Automatic ORI Detection**: Uses GC skew analysis and DnaA box detection
- **Multiple Cloning Site Construction**: Builds MCS with specified enzymes
- **Marker Integration**: Incorporates antibiotic resistance and selection markers
- **Restriction Site Management**: Removes unwanted sites
- **Broad Host Range Support**: Includes origins for diverse bacterial hosts
- **Comprehensive Testing**: Full test suite included

## Quick Start

```bash
python plasmid_designer.py sample_pUC19.fa sample_Design_pUC19.txt output.fa
```

## Installation

Requirements: Python 3.7+
```bash
# No installation needed - uses Python standard library only
```

## Usage

```bash
python plasmid_designer.py <input.fa> <design.txt> <output.fa> [markers.tab]
```

## Design File Format

```
BamHI_site, BamHI
HindIII_site, HindIII
AmpR_gene, Ampicillin
ori_pMB1, High_Copy_Replication
```

## Testing

```bash
python test_plasmid_designer.py
```

## Documentation

See full documentation in this README for:
- Detailed usage examples
- File format specifications
- Algorithm descriptions
- Troubleshooting guide
- API reference

## Based on Scientific Research

Implementation based on "Broad host range plasmids" (Jain & Srivastava, FEMS 2013)
