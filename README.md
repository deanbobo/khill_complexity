# khill_complexity.jl

A Julia implementation of **K-Hill / k-mer complexity analysis** for genome collections.

Given a directory of genomes (FASTA or gzipped FASTA), this script:

- extracts sequences
- samples / hashes k-mers
- computes per-genome and ensemble **entropy** and **effective diversity**
- computes **Kullback–Leibler divergence**–based beta-entropy and beta-diversity
- summarizes everything into `ent_hills.out` and `structure.out` for downstream plotting and analysis.

The implementation uses **MurmurHash3 (x64_128)** with a fixed seed for reproducible k-mer hashing.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/deanbobo/khill_complexity.git
cd khill_complexity
