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

## Installation & Implementation

### 1. Clone the repository

```bash
git clone https://github.com/deanbobo/khill_complexity.git
cd khill_complexity
```

### 2. Ensure you have Julia installed
```bash
julia --version
```
Julia >= 1.9 required. (tested with 1.11)

 ### 3. Install Julia package dependencies
 ```bash
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'
```

 ### 4. Put your genomes in a directory
 ```text
  input_genomes/
     sample1.fa
     sample2.fa
     sample3.fa.gz
```

**Notes**
- Files can be FASTA (.fa, .fasta) or gzipped FASTA (.fa.gz, .fasta.gz).

- Filenames (without .gz) are treated as genome IDs.

- Multi-FASTA is supported: all sequences per file are concatenated per “genome”.


### 5. Run khill_complexity.jl
```bash
JULIA_NUM_THREADS=8 julia khill_complexity.jl ...
```
**Command line arguments**
| Flag | Required | Type     | Meaning |
|------|----------|----------|---------|
| `-i` | Yes      | String   | Input directory containing genomes (`.fa`, `.fa.gz`) |
| `-o` | Yes      | String   | Output directory (created if it doesn’t exist) |
| `-k` | Yes      | String   | K-mer size or range, e.g. `"9"` or `"5-50"` |
| `-t` | Yes      | String   | Taxon identifier (label used in `structure.out`) |
| `-n` | Yes      | Int      | Sampling rate (see below) |
| `-c` | No       | Int      | Canonical kmers: `0` = as-is, `1` = canonical (min of k-mer and reverse-complement) |
| `-s` | No       | Float64  | K-mer limit per genome per k (default `1e99` = effectively no limit) |
| `-a` | No       | Int      | Presence/absence mode: `0` = count multiplicity, `1` = presence/absence per genome |
| `-x` | No       | Int      | Hashing mode: `1` = hashed (Murmur3), `2` = raw k-mers (no hashing) |



### 6. Outputs
```text
outdir/
  ent_hills.out
  structure.out
  9/
    stats.out
    sample1.betaent
    sample2.betaent
    ...
  10/
    stats.out
    ...
  ...
```


structure.out

At the top level: outdir/structure.out.

This file is a single line with 15 fields:

| Column | Name            | Meaning |
|--------|------------------|---------|
| 1  | `taxa`        | Label from `-t` |
| 2  | `totlen`      | Total genome length (sum of all bases across genomes) |
| 3  | `meanlen`     | Mean genome length |
| 4  | `stdv`        | Standard deviation of genome length |
| 5  | `deltaH_maxk` | ΔH at the largest k-mer size |
| 6  | `exent`       | Overall excess entropy across k |
| 7  | `transinfo`   | Cumulative transient information |
| 8  | `pred`        | Predictability measure (sum of Δ²H) |
| 9  | `maxh`        | Maximum ensemble entropy across k |
| 10 | `maxeffk`     | Maximum effective k-mer richness (`exp(maxh)`) |
| 11 | `maxeffg`     | Maximum effective number of genomes (`exp(betaEnt)`) |
| 12 | `minnh`       | Minimum normalized entropy across k |
| 13 | `minnhk`      | k-mer size at which `minnh` occurs |
| 14 | `mind2h`      | Minimum Δ²H across k |
| 15 | `mind2hk`     | k-mer size at which `mind2h` occurs |








