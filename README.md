# SVJedi-graph : long-read SV genotyper with a variation graph

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](https://anaconda.org/bioconda/svjedi-graph)

SVJedi-graph is a structural variation (SV) genotyper for long read data. It takes as input a variant file (VCF), a reference genome (fasta) and a long read file (fasta/fastq) and outputs the initial variant file with an additional column containing genotyping information (VCF).

SVjedi-graph is based on a representation of the genome and the different SV alleles in a variation graph. After building this variation graph from the reference genome sequence and the input variant file, long reads are mapped on this graph using minigraph[^1]. Then it estimates the genotype of each variant in a given individual sample based on allele-specific alignment counts. 

Currently, SVJedi-graph can genotype five types of SVs: deletions, insertions, duplications, inversions and translocations (intra- and inter-chromosomal).


[^1]: Li, H., Feng, X. & Chu, C. The design and construction of reference pangenome graphs with minigraph. Genome Biol 21, 265 (2020). https://doi.org/10.1186/s13059-020-02168-z

## Installation

SVJedi-graph requires :

* Python (3.8.13 or higher) 
* [minigraph](https://github.com/lh3/minigraph)

### With Conda

```bash
conda install -c bioconda svjedi-graph
```

### Or

```bash
git clone https://gitlab.inria.fr/sromain/svjedi-graph.git
```

## Usage

```bash
./svjedi-graph.py -v <inputVCF> -r <refFA> -q <longreadsFQ> [ -p <output_prefix> -t <threads> ]
``` 

### Input VCF requirements

For all variants, the `SVTYPE` tag must be present in the `INFO` field (`SVTYPE=DEL` or `SVTYPE=INS` or `SVTYPE=INV` or `SVTYPE=BND`). Insertions need to be sequence-resolved with the full inserted sequence characterized and reported in the ALT field of the VCF file. As duplications are a special case of insertions, SVJedi-graph supports also duplications, as long as their duplicated sequence is characterized and reported similarly to insertions. More details are given in [SV representation in VCF](#SV-representation-in-VCF).


### Test with a small dataset


```bash
./svjedi-graph.py -v test_data/tiny_DEL.vcf -r test_data/tiny_GRCh37_chr22.fa -q tiny_sim_reads.fastq.gz -p tiny_DEL
``` 

The genotyped VCF file obtained (`tiny_DEL_genotype.vcf`) should be identical as the one in the `test_data` folder.


### Parameters

|:---:|:---:|:---|
|`-v`|`--vcf`|VCF file containing the set of SVs to genotype.|
|`-r`| `--ref` | FASTA file containing the reference genome (on which the SVs have been identified).|
|`-q`| `--reads` | FASTQ file containing the long reads used to genotype.|
|`-p`| `--prefix` | Prefix of output files.|
|`-t`| `--threads` | Number of threads to use for the mapping step.|
|`-ms`| `--minsupport` | Minimum number of alignments to genotype a SV (default: 3>=).|

* `-v` `--vcf`  VCF file containing the set of SVs to genotype.
* `-r` `--ref`  FASTA file containing the reference genome (on which the SVs have been identified).
* `-q` `--reads`  FASTQ file containing the long reads used to genotype.
* `-p` `--prefix`  Prefix of output files.
* `-t` `--threads`  Number of threads to use for the mapping step.
* `-ms` `--minsupport`  Minimum number of alignments to genotype a SV (default: 3>=).


### Output files

Main output file:

* `<prefix>_genotype.vcf`  Genotyped SVs set in VCF format.

Intermediate output files:

* `<prefix>.gfa`           Variation graph in [GFA format](https://github.com/GFA-spec/GFA-spec).
* `<prefix>.gaf`           Mapping results from minigraph in [GAF format](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf).
* `<prefix>_informative_aln.json`   Json dictionnary of read supports for each input SV's alleles.



### SV representation in VCF

Here are the information needed for SVJedi-graph to genotype the following SV types. All variants must have the ```CHROM``` and ```POS``` fields defined, with the chromosome names in the reference genome file and variant file that must be the same. The `SVTYPE` tag must be present in the INFO field (`SVTYPE=DEL` or `SVTYPE=INS` or `SVTYPE=INV` or `SVTYPE=BND`). Then additional information is required according to SV type:

- Deletion
	- ```INFO``` field must contain ```SVTYPE=DEL```
	- ```INFO``` field must contain ```END=pos``` (with `pos` being the end position of the deleted segment)
	
- Insertion
	- ```INFO``` field must contain ```SVTYPE=INS```
	- ```ALT``` field must contain the sequence of the insertion 

- Duplication
	- must be defined as an insertion event whith `CHR` and `POS` corresponding to the position of insertion of the novel copy
	- ```INFO``` field must contain ```SVTYPE=INS```
	- ```ALT``` field must contain the sequence of the duplication 
	
- Inversion
	- ```INFO``` field must contain ```SVTYPE=INV```
	- ```INFO``` field must contain ```END=pos``` tag, with `pos` being the second breakpoint position

- Intra-chromosomal translocation
	- ```INFO``` field must contain ```SVTYPE=BND```
	- ```ALT``` field must be formated as: ```t[pos[```, ```t]pos]```, ```]pos]t``` or ```[pos[t```, with `pos` indicating the second breakpoint position and brackets directions indicating which parts of the two chromosomes should be joined together 

## Citation

Sandra Romain, Claire Lemaitre, SVJedi-graph: improving the genotyping of close and overlapping structural variants with long reads using a variation graph, Bioinformatics, Volume 39, Issue Supplement_1, June 2023, Pages i270–i278, https://doi.org/10.1093/bioinformatics/btad237

## Contact

SVJedi-graph is a [Genscale](https://team.inria.fr/genscale/) tool developed by Sandra Romain and Claire Lemaitre. For any bug report or feedback, please use the Github Issues form.

---
