import pytest, numpy, copy

from pathlib import Path

import gumpy

import piezo

TEST_CASE_DIR = "tests/test-cases/"

reference_genome=gumpy.Genome(genbank_file="config/NC_004148.2.gbk",name="NC_004148.2")

resistance_catalogue=piezo.ResistanceCatalogue( input_file="config/NC_004148.2_RSU_catalogue.csv",
                                            gumpy_genome=reference_genome,
                                            catalogue_name="TEST" )


print(resistance_catalogue.gene_list)

def test_catalogue__init__():

    assert resistance_catalogue.catalogue_name=="TEST"

    assert resistance_catalogue.number_rows==12

    assert resistance_catalogue.drug_list==['DRUG_A','DRUG_B']

    assert resistance_catalogue.gene_list==['M2']

def test_catalogue_prediction_snps():

    # check a row in the catalogue
    assert resistance_catalogue.predict("M2_L73P")=={'DRUG_A': 'R', 'DRUG_B': 'U'}

    # check a synonymous mutation has no effect
    assert resistance_catalogue.predict("M2_L73L")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    # check hitting a wildtype row
    assert resistance_catalogue.predict("M2_G74I")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert resistance_catalogue.predict("M2_t-15c")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # bad prediction
    with pytest.raises(Exception):
        assert resistance_catalogue.predict("M2_L73P")=={'DRUG_A': 'R', 'DRUG_B': 'R'}

    # incorrect amino acid in the ref position
    with pytest.raises(Exception):
        assert resistance_catalogue.predict("M2_K73P")

    # incorrect amino acid in the alt position
    with pytest.raises(Exception):
        assert resistance_catalogue.predict("M2_L73B")

    # gene not present in genome
    with pytest.raises(Exception):
        assert resistance_catalogue.predict("M3_K73P")

    # gene not present in genome
    with pytest.raises(Exception):
        assert resistance_catalogue.predict("M3_K73P")

    # badly formed gene_mutation
    with pytest.raises(Exception):
        assert resistance_catalogue.predict("M3_K73P_3")
        assert resistance_catalogue.predict("M3-K73P")

def test_catalogue_prediction_indels():

    # check a general indel
    assert resistance_catalogue.predict("M2_300_indel")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check an insertion
    assert resistance_catalogue.predict("M2_300_ins")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert resistance_catalogue.predict("M2_-10_ins")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a deletion
    assert resistance_catalogue.predict("M2_300_del")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check numbered insertions and deletions
    assert resistance_catalogue.predict("M2_300_ins_5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert resistance_catalogue.predict("M2_300_del_5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert resistance_catalogue.predict("M2_300_indel_5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert resistance_catalogue.predict("M2_300_indel_-5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a specific insertion
    assert resistance_catalogue.predict("M2_300_ins_acct")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check these hit the frameshift rule
    assert resistance_catalogue.predict("M2_300_ins_6")=={'DRUG_A': 'U', 'DRUG_B': 'R'}
    assert resistance_catalogue.predict("M2_300_del_6")=={'DRUG_A': 'U', 'DRUG_B': 'R'}
    assert resistance_catalogue.predict("M2_300_ins_acc")=={'DRUG_A': 'U', 'DRUG_B': 'R'}

    # badly formed INDELs
    with pytest.raises(Exception):
        assert resistance_catalogue.predict("M2_300_insdel")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
