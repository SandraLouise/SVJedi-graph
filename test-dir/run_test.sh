#! /bin/bash

# look for svjedi-graph.py (in PATH for conda install, in parent directory for normal install)
if [ -x "$(command -v svjedi-graph.py)" ] 
then 
bin="svedi-graph.py"
elif [ -f "../svjedi-graph.py" ] 
then
bin="../svjedi-graph.py"
else
echo "could not find svjedi-graph.py"
exit 1
fi


RETVAL=0


Nbthreads=1
REF="reference_genome.fasta"
FASTQ="simulated_reads.fastq.gz"

VCF="test.vcf"
pref="test"

svjedi-graph.py -v $VCF -r $REF -q $FASTQ -p $pref -t $Nbthreads  1> /dev/null 2>&1

outfile=$pref"_genotype.vcf"

# Comparing genotypes

./contingency_table.py $VCF $outfile > ${outfile}".eval"

diff ${outfile}".eval" expected_genotype.vcf.eval 1> /dev/null 2>&1
var=$?

if [ $var -eq 0 ]
then
echo "SVJedi-graph test : PASS"
else
echo "SVJedi-graph test : FAILED"
RETVAL=1
fi


# Looking in more details, are VCF lines identical (and in the same order) ?

diff --ignore-matching-lines="^#"  ${outfile} expected_genotype.vcf 1> /dev/null 2>&1
var=$?

echo "-----------------"
echo "Details:"
if [ $var -eq 0 ]
then
echo "VCF lines are identical"
else
echo "Genotypes are correct but VCF lines differ"
fi




# Cleaning
rm -f ${pref}.gaf
rm -f ${pref}.gfa
rm -f ${pref}_informative_aln.json
rm -f svs_edges.json
rm -f ignored_svs.vcf

rm -f ${pref}_genotype.vcf
rm -f ${pref}_genotype.vcf.eval


# for Jenkins CI platform, we need an exit code: PASS (0) vs. FAILED (1)
exit $RETVAL
