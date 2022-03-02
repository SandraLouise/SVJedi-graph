# SVJedi-graph : long-read SV genotyper with variant graph

SVJedi-graph constructs a variation graph to represent a set of SVs from a VCF file, which is then used as a base to map long reads on with GraphAligner [^1].

[^1]: Rautiainen, M., Marschall, T. GraphAligner: rapid and versatile sequence-to-graph alignment. Genome Biol 21, 253 (2020). https://doi.org/10.1186/s13059-020-02157-2.

## Installation

SVJedi-graph requires only python and GraphAligner to run. Currently works best with GraphAligner v1.0.12 (https://github.com/maickrau/GraphAligner/releases/tag/v1.0.12). Simplest way to install GraphAligner is to use a conda environment.

```bash
conda install -c bioconda graphaligner=1.0.12
```

## Run

Make sure to activate your environment if you installed GraphAligner _via_ Conda.

```bash
python3 svjedi-graph.py -v <inputVCF> -r <refFA> -q <longreadsFQ> -p <output_prefix> -t <threads>
```


### Parameters


* `-v`  VCF file containing the set of SVs to genotype.
* `-r`  FASTA file containing the reference genome (on which the SVs have been identified).
* `-q`  FASTQ file containing the long reads used to genotype.
* `-p`  Prefix of output files.
* `-t`  Number of threads to use for the mapping step.

### Output files

* `<prefix>.gfa`           Variation graph in [GFA format](https://github.com/GFA-spec/GFA-spec).
* `<prefix>.gaf`           Mapping results from GraphAligner in [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf).
* `<prefix>_genotype.vcf`  Genotyped SVs set in VCF format.

