# SVJedi-graph : long-read SV genotyper with variant graph

SVJedi-graph constructs a variation graph to represent a set of SVs from a VCF file, which is then used as a base to map long reads on with GraphAligner [^1]. Currently works best with GraphAligner v1.0.12 (https://github.com/maickrau/GraphAligner/releases/tag/v1.0.12).

[^1]: Rautiainen, M., Marschall, T. GraphAligner: rapid and versatile sequence-to-graph alignment. Genome Biol 21, 253 (2020). https://doi.org/10.1186/s13059-020-02157-2.

## Installation

SVJedi-graph requires only python and GraphAligner to run.

```bash
...
```

## Usage

```bash
...
```

Parameters:

* `-v`  VCF file containing the set of SVs to genotype.
* `-r`  FASTA file containing the sequence of the reference genome (on which the SVs have been identified).
* `-o`  Prefix of output files.

## Output files

* `<prefix>_genotyped.vcf`  Genotyped SVs set in VCF format.
* `<prefix>.gfa`            Variation graph in GFA format.
