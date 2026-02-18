# Plasmid Designer Tool

A bioinformatics tool for designing plasmids from genomic sequences. Automatically detects origins of replication (ORI) and builds functional plasmids with specified restriction sites and selectable markers.

## Requirements

Python 3.7+ — no additional installation needed.

## Quick Start

```bash
python3 plasmid_designer.py pUC19.fa Design_pUC19.txt Output.fa markers.tab
```

## File Formats

**Input FASTA** (`pUC19.fa`):
```
>Header
ATGCATGCATGC...
```

**Design file** (`Design_pUC19.txt`) — one component per line, `Name, Marker`:
```
BamHI_site, BamHI
HindIII_site, HindIII
AmpR_gene, Ampicillin
ori_pMB1, High_Copy_Replication
```

**Markers database** (`markers.tab`) — tab-separated:
```
Name    Type    Sequence    Description
```

## Usage

```bash
# Design a plasmid
python3 plasmid_designer.py <input.fa> <design.txt> <output.fa> [markers.tab]

# Run tests
python3 test_plasmid_designer.py

# Analyze output
python3 plasmid_analyzer.py <output.fa>
```

## Features

- **ORI Detection** — GC skew analysis, DnaA box pattern matching, AT content
- **MCS Construction** — builds multiple cloning site from specified enzymes
- **Marker Integration** — antibiotic resistance and selection markers
- **Restriction Site Management** — automatically removes sites not in design file
- **Broad Host Range** — supports diverse bacterial hosts

## Common Tasks

| Goal | How |
|------|-----|
| Add restriction sites | Add `SiteName, EnzymeName` to design file |
| Add antibiotic resistance | Add `AmpR_gene, Ampicillin` to design file |
| Remove a site | Simply omit it from the design file |
| Verify site removal | `python3 plasmid_analyzer.py Output.fa \| grep EcoRI` |

## Project Files

```
plasmid_designer.py       # Main tool
test_plasmid_designer.py  # Test suite
plasmid_analyzer.py       # Output analyzer
markers.tab               # Marker database
pUC19.fa                  # Example input
Design_pUC19.txt          # Example design
```

## Based On

Jain & Srivastava, "Broad host range plasmids", *FEMS Microbiology Letters*, 2013.
