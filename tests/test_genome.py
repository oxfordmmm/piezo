import pytest, numpy, copy

from pathlib import Path

from piezo import Genome

TEST_CASE_DIR = "tests/test-cases/"

reference=Genome(genbank_file="config/H37rV_v3.gbk")

def test_Genome_instantiate_genbank():

    # check that the M. tuberculosis H37rV genome is the right length
    assert reference.length==4411532

    # check the species is stored correctly
    assert reference.organism=='Mycobacterium tuberculosis H37Rv'

    # check that the sequence starts and ends as we expect
    assert reference.sequence[0]=='t'
    assert reference.sequence[-1]=='g'

    assert len(reference.gene_names)==3909

    # try a normal gene
    assert reference.gene_is_reverse['rpoB'] is False
    assert reference.gene_type['rpoB']=="GENE"
    mask=reference.gene=="rpoB"
    sequence=reference.sequence[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="ttg"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="taa"

    # try a reverse complement gene
    assert reference.gene_is_reverse['katG'] is True
    assert reference.gene_type['katG']=="GENE"
    mask=reference.gene=="katG"
    sequence=reference.sequence[mask]
    first_codon="".join(i for i in sequence[:3])
    assert first_codon=="tca"
    last_codon="".join(i for i in sequence[-3:])
    assert last_codon=="cac"


reference2=Genome(fasta_file="config/H37rV_v3.fasta.gz")

def test_Genome_instantiate_fasta():

    # check that the M. tuberculosis H37rV genome is the right length
    assert reference2.length==4411532

    # check the species is stored correctly
    assert reference2.organism=='Mycobacterium tuberculosis H37Rv'

    # check that the sequence starts and ends as we expect
    assert reference2.sequence[0]=='t'
    assert reference2.sequence[-1]=='g'


def test_Genome_gbk_fasta_identical():

    assert reference.length==reference2.length

    assert numpy.array_equal(reference.sequence,reference2.sequence)

    assert numpy.array_equal(reference.diploid,reference2.diploid)


def test_Genome___repr__():

    assert reference.__repr__()=='NC_000962.3\nMycobacterium tuberculosis H37Rv\nReference Genome\n4411532 bases\nttg...tcg'


def test_Genome___sub__():

    sample=copy.deepcopy(reference)

    sample.sequence[2]='t' # remember that the genome is 1-based, but the numpy array is 0-based

    (original_bases,indices,new_bases)=reference-sample

    assert original_bases[0]=='g'
    assert new_bases[0]=='t'
    assert indices[0]==3

def test_Genome_calculate_snp_distance():

    sample=copy.deepcopy(reference)

    sample.sequence[2]='t' # remember that the genome is 1-based, but the numpy array is 0-based
    assert sample.calculate_snp_distance(reference)==1

    # reverse the change
    sample.sequence[2]='g'
    assert sample.calculate_snp_distance(reference)==0

    sample.sequence[2]='t'
    sample.sequence[3]='t'
    assert sample.calculate_snp_distance(reference)==2

def test_Genome_apply_vcf():

    sample_03=copy.deepcopy(reference)
    sample_03.apply_vcf_file(vcf_file=TEST_CASE_DIR+"03.vcf",ignore_status=True,ignore_filter=True)
    (original_bases,indices,new_bases)=reference-sample_03
    assert original_bases[0]=='c'
    assert indices[0]==761155
    assert new_bases[0]=='t'

    sample_04=copy.deepcopy(reference)
    sample_04.apply_vcf_file(vcf_file=TEST_CASE_DIR+"04.vcf",ignore_status=True,ignore_filter=True)
    (original_bases,indices,new_bases)=reference-sample_04
    assert original_bases[0]=='c'
    assert indices[0]==2155168
    assert new_bases[0]=='g'
