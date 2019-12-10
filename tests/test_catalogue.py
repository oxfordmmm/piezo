import pytest

import piezo

test=piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv")

def test_catalogue__init__():

    assert test.catalogue.genbank_reference=="NC_004148.2"

    assert test.catalogue.name=="TEST"

    assert test.catalogue.version=="v1.0"

    assert test.catalogue.values==["R","F","U","S"]

    assert test.catalogue.grammar=="GARC1"

    assert test.catalogue.number_rows==20

    assert test.catalogue.drugs==['DRUG_A','DRUG_B']

    assert test.catalogue.genes==['M2']

    assert test.catalogue.gene_lookup == {'M2':['DRUG_A','DRUG_B']}

    assert test.catalogue.drug_lookup == {'DRUG_A':['M2'],'DRUG_B':['M2']}
#
def test_catalogue_prediction_snps():

    # check a row in the catalogue
    assert test.predict("M2_L73P")=={'DRUG_A': 'R', 'DRUG_B': 'U'}

    # check a synonymous mutation has no effect
    assert test.predict("M2_L73L")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    # check a row in the catalogue
    assert test.predict("M2_L73Z")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    # check hitting a wildtype row
    assert test.predict("M2_G74I")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert test.predict("M2_G74Z")=={'DRUG_A': 'R', 'DRUG_B': 'S'}

    assert test.predict("M2_G74P")=={'DRUG_A': 'R', 'DRUG_B': 'U'}

    assert test.predict("M2_t-15c")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert test.predict("M2_t-15z")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    # check that a gene not in the catalogue simply returns an "S"
    assert test.predict("N1_S2T")=="S"

    # bad prediction
    with pytest.raises(Exception):
        assert test.predict("M2_L73P")=={'DRUG_A': 'R', 'DRUG_B': 'R'}

    # incorrect amino acid in the alt position
    with pytest.raises(Exception):
        assert test.predict("M2_L73B")

    # badly formed gene_mutation
    with pytest.raises(Exception):
        assert test.predict("M3_K73P_3")
        assert test.predict("M3-K73P")
        assert test.predict("M2_K73t")
        assert test.predict("M2_a-10a")
        assert test.predict("M2_a-10A")
        assert test.predict("M2_T-10a")
#
def test_catalogue_prediction_indels():

    # check a general indel
    assert test.predict("M2_300_indel")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check an insertion
    assert test.predict("M2_300_ins")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert test.predict("M2_-10_ins")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a deletion
    assert test.predict("M2_300_del")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check numbered insertions and deletions
    assert test.predict("M2_300_ins_5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert test.predict("M2_300_del_5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a specific insertion
    assert test.predict("M2_300_ins_acct")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert test.predict("M2_100_indel")=={'DRUG_A': 'R', 'DRUG_B': 'U'}
    assert test.predict("M2_100_ins")=={'DRUG_A': 'R', 'DRUG_B': 'U'}
    assert test.predict("M2_100_ins_2")=={'DRUG_A': 'R', 'DRUG_B': 'U'}
    assert test.predict("M2_100_ins_ac")=={'DRUG_A': 'R', 'DRUG_B': 'U'}


    # check these hit the frameshift rule
    assert test.predict("M2_100_ins_3")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test.predict("M2_100_ins_act")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test.predict("M2_300_ins_6")=={'DRUG_A': 'U', 'DRUG_B': 'R'}
    assert test.predict("M2_300_del_6")=={'DRUG_A': 'U', 'DRUG_B': 'R'}
    assert test.predict("M2_300_ins_acc")=={'DRUG_A': 'U', 'DRUG_B': 'R'}

    # badly formed INDELs
    with pytest.raises(Exception):
        assert test.predict("M2_300_indel_5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
        assert test.predict("M2_300_indel_-5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
        assert test.predict("M2_300_insdel")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
        assert test.predict("M2_300_del_act")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
        assert test.predict("M2_300_del_ACT")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
