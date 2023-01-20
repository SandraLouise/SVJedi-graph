# SVJedi-graph : long-read SV genotyper with variant graph

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/svjedi-graph)

SVJedi-graph constructs a variation graph to represent a set of SVs from a VCF file, which is then used as a base to map long reads on with minigraph[^1].

[^1]: Li, H., Feng, X. & Chu, C. The design and construction of reference pangenome graphs with minigraph. Genome Biol 21, 265 (2020). https://doi.org/10.1186/s13059-020-02168-z

## Installation

SVJedi-graph requires only python (3.8.13 or higher) and minigraph to run.

### With Conda

```bash
conda install -c bioconda svjedi-graph
```

### Or

```bash
git clone https://gitlab.inria.fr/sromain/svjedi-graph.git
```

## Run

```bash
python3 svjedi-graph.py -v <inputVCF> -r <refFA> -q <longreadsFQ> [ -p <output_prefix> -t <threads> ]
``` 

### Parameters

* `-v`  VCF file containing the set of SVs to genotype.
* `-r`  FASTA file containing the reference genome (on which the SVs have been identified).
* `-q`  FASTQ file containing the long reads used to genotype.
* `-p`  Prefix of output files.
* `-t`  Number of threads to use for the mapping step.

### Output files

* `<prefix>.gfa`           Variation graph in [GFA format](https://github.com/GFA-spec/GFA-spec).
* `<prefix>.gaf`           Mapping results from minigraph in [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf).
* `<prefix>_informative_aln.json`   Json dictionnary of read supports for each input SV's alleles.
* `<prefix>_genotype.vcf`  Genotyped SVs set in VCF format.

## Contact

SVJedi is a [Genscale](https://team.inria.fr/genscale/) tool developed by Sandra Romain and Claire Lemaitre. For any bug report or feedback, please use the Issues form.

---
