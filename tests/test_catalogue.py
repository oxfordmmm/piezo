import pytest

import piezo

test=piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv")

def test_catalogue__init__():

    assert test.catalogue.genbank_reference=="NC_004148.2"

    assert test.catalogue.name=="TEST"

    assert test.catalogue.version=="v1.0"

    assert test.catalogue.values==["R","F","U","S"]

    assert test.catalogue.grammar=="GARC1"

    assert test.catalogue.number_rows==36

    assert test.catalogue.drugs==['DRUG_A','DRUG_B']

    assert test.catalogue.genes==['M2']

    assert test.catalogue.gene_lookup == {'M2':['DRUG_A','DRUG_B']}

    assert test.catalogue.drug_lookup == {'DRUG_A':['M2'],'DRUG_B':['M2']}
#
def test_catalogue_prediction_snps():

    # check a row in the catalogue
    assert test.predict("M2@L73P")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a synonymous mutation has no effect
    assert test.predict("M2@L73L")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    # check a het
    assert test.predict("M2@L73Z")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    # check a filter fail that will hit the default S rule
    assert test.predict("M2@L73O")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    # check a null that will hit the default S rule
    assert test.predict("M2@L73X")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    # check a filter fail that hits a specific rule for DRUG_A
    assert test.predict("M2@G74O")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    # check a null that hits a specific rule for DRUG_A
    assert test.predict("M2@G74X")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    # check hitting a wildtype row
    assert test.predict("M2@G74I")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert test.predict("M2@G74Z")=={'DRUG_A': 'R', 'DRUG_B': 'S'}

    assert test.predict("M2@G74P")=={'DRUG_A': 'R', 'DRUG_B': 'U'}

    assert test.predict("M2@G74!")=={'DRUG_A': 'S', 'DRUG_B': 'U'}

    assert test.predict("M2@t-15c")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert test.predict("M2@t-15z")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    assert test.predict("M2@t-15o")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    assert test.predict("M2@t-15x")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    # check that a gene not in the catalogue simply returns an "S"
    assert test.predict("N1@S2T")=="S"

    # bad prediction
    with pytest.raises(Exception):
        assert test.predict("M2@L73P")=={'DRUG_A': 'R', 'DRUG_B': 'R'}

    # incorrect amino acid in the alt position
    with pytest.raises(Exception):
        assert test.predict("M2@L73B")

    # badly formed gene_mutation
    with pytest.raises(Exception):
        assert test.predict("M3@K73P_3")
        assert test.predict("M3@K73P")
        assert test.predict("M2@K73t")
        assert test.predict("M2@a-10a")
        assert test.predict("M2@a-10A")
        assert test.predict("M2@T-10a")
#
def test_catalogue_prediction_indels():

    # check a general indel
    assert test.predict("M2@300_indel")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check an insertion
    assert test.predict("M2@300_ins")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert test.predict("M2@-10_ins")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a deletion
    assert test.predict("M2@300_del")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check numbered insertions and deletions
    assert test.predict("M2@300_ins_6")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert test.predict("M2@300_del_6")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a specific insertion
    assert test.predict("M2@300_ins_acc")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert test.predict("M2@100_indel")=={'DRUG_A': 'R', 'DRUG_B': 'U'}
    assert test.predict("M2@100_ins")=={'DRUG_A': 'R', 'DRUG_B': 'U'}
    assert test.predict("M2@100_ins_3")=={'DRUG_A': 'R', 'DRUG_B': 'U'}
    assert test.predict("M2@100_ins_act")=={'DRUG_A': 'R', 'DRUG_B': 'U'}


    # check these hit the frameshift rule
    assert test.predict("M2@100_ins_1")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test.predict("M2@100_ins_acta")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test.predict("M2@100_ins_400")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test.predict("M2@100_ins_actgactg")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test.predict("M2@300_ins_61")=={'DRUG_A': 'U', 'DRUG_B': 'R'}
    assert test.predict("M2@300_del_65")=={'DRUG_A': 'U', 'DRUG_B': 'R'}
    assert test.predict("M2@300_ins_acctt")=={'DRUG_A': 'U', 'DRUG_B': 'R'}

    # badly formed INDELs
    with pytest.raises(Exception):
        assert test.predict("M2@300_indel_5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
        assert test.predict("M2@300_indel_-5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
        assert test.predict("M2@300_insdel")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
        assert test.predict("M2@300_del_act")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
        assert test.predict("M2@300_del_ACT")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
