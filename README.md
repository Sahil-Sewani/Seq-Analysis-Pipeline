# Seq-Analysis-Pipeline
## Program Summary and Interconnection

## 1. create_paired_aaseqs.py
### Overview
This script takes a file containing pairs of eukaryotic species and aligns their DNA sequences based on aligned amino acid sequences. It generates new aligned DNA sequences.

### Interconnection
- Utilizes functions from `AAtoDNA.py` for aligning DNA sequences based on protein alignment.
- Functions in `AAtoDNA.py` are designed for general-purpose amino acid to DNA alignment.

## 2. align_aa_seqs.py
### Overview
This program reads pairs of aligned DNA sequences and calculates various values, outputting results to a file. It considers different sets of degenerate sites, transitions, transversions, and calculates Kappa values.

### Interconnection
- Uses functions from `calculate_from_pairs.py` for filtering, counting, and calculating values.
- Shared functions are designed for general-purpose sequence analysis.

## 3. AAtoDNA.py
### Overview
This script converts aligned amino acid sequences to aligned DNA sequences by inserting corresponding codons. It is used by `create_paired_aaseqs.py` and `align_aa_seqs.py` for aligning DNA sequences based on protein alignment.

### Interconnection
- Invoked by `create_paired_aaseqs.py` and `align_aa_seqs.py` for the core functionality of aligning DNA sequences based on amino acid alignment.

## 4. calculate_from_pairs.py
### Overview
This program reads pairs of aligned DNA sequences and calculates various values, outputting results to a file. It shares some functions with `align_aa_seqs.py` for filtering, counting, and calculating values.

### Interconnection
- Functions are shared with `align_aa_seqs.py` for common tasks such as filtering and calculating various counts.
- Both scripts perform calculations on the same input files but focus on different aspects of the analysis.

## Execution
- All programs can be executed from the command line, taking the species pair file as input.
- The interconnection allows them to leverage each other's functionalities for a comprehensive analysis of aligned DNA sequences.

