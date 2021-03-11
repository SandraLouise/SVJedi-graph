# SVJedi-graph : long-red SV genotyper with variant graph

Graph construction with VG toolkit (v?).
Long-read mapping on graph with GraphAligner (v1.0.12).

### Current workflow

**1. Variant data prep**

```
python3 extract_variant_data.py -v <vcf_input> -r <genome_ref_input> 
```

Output : reference sequence for all SV regions in vcf input (`reference_subsequences.fa`), converted SV positions on their region's reference sequence (`converted_variants.vcf`).


**2. Graph construction**

```
vg construct -r reference_subsequences.fa -v converted_variants.vcf -m 5000 -S -a > variant_graph.vg
vg view --gfa variant_graph.vg > variant_graph.gfa
```

Output : variant graph in vg (`variant_graph.vg`) and gfa format (`variant_graph.gfa`).

**3. Long-read mapping**

```
GraphAligner -g variant_graph.vg -f <long_read_input> -a aln.gaf -x vg > log.txt
```

Output : alignment results in gaf format (`aln.gaf`), log of GraphAligner run (`log.txt`).

**4. Alignment filtering**

```
python3 filter_aln.py -a aln.gaf -g variant_graph.gfa
```

Output : dictionary of informative alignments on SVs in json format (`informative_aln.json`).

**5. Genotype prediction**

```
python3 svjedi_genotype_prediction.py -d informative_aln.json -v <vcf_input>
```

Output : vcf file with predicted genotype for reads (`genotype_results.txt`).

Optional arguments: 

* `-o` for output file name
