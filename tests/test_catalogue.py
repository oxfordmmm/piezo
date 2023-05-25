from collections import defaultdict
import pytest

import piezo


@pytest.mark.parametrize(
                        "test_catalogue,genes,gene_lookup", 
                        [
                            (
                                piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"), 
                                ['M2', 'N', 'MULTI'],
                                {'M2':['DRUG_A','DRUG_B'], 'MULTI': ['DRUG_A', 'DRUG_B'], 'N': ['DRUG_A']}
                            ), 
                            (
                                piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv", prediction_subset_only=True), 
                                ['M2', 'MULTI'],
                                {'M2':['DRUG_A','DRUG_B'], 'MULTI': ['DRUG_A', 'DRUG_B']}
                            )
                        ]
                        )
def test_catalogue__init__(test_catalogue, genes, gene_lookup):

    assert test_catalogue.catalogue.genbank_reference=="NC_004148.2"

    assert test_catalogue.catalogue.name=="TEST"

    assert test_catalogue.catalogue.version=="v1.0"

    assert test_catalogue.catalogue.values==["R","F","U","S"]

    assert test_catalogue.catalogue.grammar=="GARC1"

    assert test_catalogue.catalogue.number_rows == 45

    assert test_catalogue.catalogue.drugs==['DRUG_A','DRUG_B']

    #Ordering of lists is annoying here, so sort both
    genes = sorted(genes)
    assert sorted(test_catalogue.catalogue.genes) == genes

    assert test_catalogue.catalogue.gene_lookup == gene_lookup

    #Convert the gene_lookup to drug_lookup
    lookup = defaultdict(list)
    for gene in gene_lookup.keys():
        for drug in gene_lookup[gene]:
            lookup[drug].append(gene)
    lookup = {drug: sorted(lookup[drug]) for drug in lookup.keys()}

    actual_lookup = {drug: sorted(test_catalogue.catalogue.drug_lookup[drug]) for drug in test_catalogue.catalogue.drug_lookup.keys()}
    assert actual_lookup == lookup

@pytest.mark.parametrize("test_catalogue", [piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"), piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv", prediction_subset_only=True)])
def test_catalogue_prediction_snps(test_catalogue):

    # check a row in the catalogue
    assert test_catalogue.predict("M2@L73P")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a synonymous mutation has no effect
    assert test_catalogue.predict("M2@L73L")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    # check a het
    assert test_catalogue.predict("M2@L73Z")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    # check a filter fail that will hit the default S rule
    assert test_catalogue.predict("M2@L73O")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    # check a null that will hit the default S rule
    assert test_catalogue.predict("M2@L73X")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    # check a filter fail that hits a specific rule for DRUG_A
    assert test_catalogue.predict("M2@G74O")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    # check a null that hits a specific rule for DRUG_A
    assert test_catalogue.predict("M2@G74X")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    # check hitting a wildtype row
    assert test_catalogue.predict("M2@G74I")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert test_catalogue.predict("M2@G74Z")=={'DRUG_A': 'R', 'DRUG_B': 'S'}

    assert test_catalogue.predict("M2@G74P")=={'DRUG_A': 'R', 'DRUG_B': 'U'}

    assert test_catalogue.predict("M2@G74!")=={'DRUG_A': 'S', 'DRUG_B': 'U'}

    assert test_catalogue.predict("M2@t-15c")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert test_catalogue.predict("M2@t-15z")=={'DRUG_A': 'S', 'DRUG_B': 'S'}

    assert test_catalogue.predict("M2@t-15o")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    assert test_catalogue.predict("M2@t-15x")=={'DRUG_A': 'F', 'DRUG_B': 'S'}

    # check that a gene not in the catalogue simply returns an "S"
    assert test_catalogue.predict("N1@S2T")=="S"

    #Checking large deletions
    assert test_catalogue.predict("M2@del_1.0") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.9") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.85") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.8") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.7") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.5") == {"DRUG_A": "U", "DRUG_B": "U"}

    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@del_0.5s") == {"DRUG_A": "U", "DRUG_B": "U"}

    # bad prediction
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@L73P")=={'DRUG_A': 'R', 'DRUG_B': 'R'}

    # incorrect amino acid in the alt position
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@L73B")

    # badly formed gene_mutation
    with pytest.raises(Exception):
        assert test_catalogue.predict("M3@K73P_3")
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@K73t")
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@a-10a")
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@a-10A")
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@T-10a")

@pytest.mark.parametrize("test_catalogue", [piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"), piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv", prediction_subset_only=True)])
def test_catalogue_prediction_indels(test_catalogue):

    # check a general indel
    assert test_catalogue.predict("M2@300_indel")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check an insertion
    assert test_catalogue.predict("M2@300_ins")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert test_catalogue.predict("M2@-10_ins")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert test_catalogue.predict("M2@-10_ins_2")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a deletion
    assert test_catalogue.predict("M2@300_del")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check numbered insertions and deletions
    assert test_catalogue.predict("M2@300_ins_6")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    assert test_catalogue.predict("M2@300_del_6")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    # check a specific insertion
    assert test_catalogue.predict("M2@300_ins_acc")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

    assert test_catalogue.predict("M2@100_indel")=={'DRUG_A': 'R', 'DRUG_B': 'U'}
    assert test_catalogue.predict("M2@100_ins")=={'DRUG_A': 'R', 'DRUG_B': 'U'}
    assert test_catalogue.predict("M2@100_ins_3")=={'DRUG_A': 'R', 'DRUG_B': 'U'}
    assert test_catalogue.predict("M2@100_ins_act")=={'DRUG_A': 'R', 'DRUG_B': 'U'}


    # check these hit the frameshift rule
    assert test_catalogue.predict("M2@100_ins_1")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test_catalogue.predict("M2@100_ins_acta")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test_catalogue.predict("M2@100_ins_400")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test_catalogue.predict("M2@100_ins_actgactg")=={'DRUG_A': 'R', 'DRUG_B': 'R'}
    assert test_catalogue.predict("M2@300_ins_61")=={'DRUG_A': 'U', 'DRUG_B': 'R'}
    assert test_catalogue.predict("M2@300_del_65")=={'DRUG_A': 'U', 'DRUG_B': 'R'}
    assert test_catalogue.predict("M2@300_ins_acctt")=={'DRUG_A': 'U', 'DRUG_B': 'R'}

    # badly formed INDELs
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@300_indel_5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@300_indel_-5")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@300_insdel")=={'DRUG_A': 'U', 'DRUG_B': 'U'}
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@300_del_ACT")=={'DRUG_A': 'U', 'DRUG_B': 'U'}

@pytest.mark.parametrize("test_catalogue", [piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"), piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv", prediction_subset_only=True)])
def test_multi(test_catalogue):
    #Exact match with 2 lines of the catalogue
    assert test_catalogue.predict("M2@G74!&M2@G74X") == {"DRUG_A": "U", "DRUG_B": "R"}

    #Should hit an 'S' (and 'U' default) and an 'R'
    assert test_catalogue.predict("M2@F75T&M2@G74P") == {"DRUG_A": 'R', "DRUG_B": 'U'}

    #Should hit an 'F' and an 'R' (and 'U' default)
    assert test_catalogue.predict("M2@L73Z&M2@G74P") == {"DRUG_A": 'R', 'DRUG_B': 'U'}

    #Shouldn't hit anything other than 'S'
    assert test_catalogue.predict("M2@A5A&M2@K6K") == 'S'

@pytest.mark.parametrize("test_catalogue", [piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"), piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv", prediction_subset_only=True)])
def test_minor_population(test_catalogue):
    #COV
    #Exact match
    assert test_catalogue.predict("M2@F75V:2") == {"DRUG_A": "R", "DRUG_B": "U"}

    #Greater than
    assert test_catalogue.predict("M2@F75V:3") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@F75V:4") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@F75V:350") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@F75V:123456789123456789") == {"DRUG_A": "R", "DRUG_B": "U"}

    #Less than (so should just hit a default for now)
    assert test_catalogue.predict("M2@F75V:1") == {"DRUG_A": "U", "DRUG_B": "U"}

    #Should fail
    with pytest.raises(Exception):
        test_catalogue.predict("M2@F75V:0")
    with pytest.raises(Exception):
        test_catalogue.predict("M2@F75V:-1")

    #Multi
    assert test_catalogue.predict("M2@45_del_aaa:6&M2@F75V:5") == {"DRUG_A": 'R', 'DRUG_B': 'U'}
    assert test_catalogue.predict("M2@45_del_a:2&M2@F75V:3") == {"DRUG_A": 'R', 'DRUG_B': 'R'}

    #FRS
    test_catalogue2 = piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS-FRS.csv")

    #Exact match
    assert test_catalogue2.predict("M2@F75V:0.03") == {"DRUG_A": "R", "DRUG_B": "U"}

    #Greater than
    assert test_catalogue2.predict("M2@F75V:0.04") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue2.predict("M2@F75V:0.05") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue2.predict("M2@F75V:0.9") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue2.predict("M2@F75V:0.123456789123456789") == {"DRUG_A": "R", "DRUG_B": "U"}

    #Less than (so should just hit a default for now)
    assert test_catalogue2.predict("M2@F75V:0.01") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue2.predict("M2@F75V:0.02") == {"DRUG_A": "U", "DRUG_B": "U"}

    #Should fail
    with pytest.raises(Exception):
        test_catalogue2.predict("M2@F75V:0.0")
    with pytest.raises(Exception):
        test_catalogue2.predict("M2@F75V:-0.1")
    
    assert test_catalogue2.predict("M2@45_del_aaa:0.01&M2@F75V:0.5") == {"DRUG_A": 'R', 'DRUG_B': 'U'}
    assert test_catalogue2.predict("M2@45_del_a:0.12&M2@F75V:0.03") == {"DRUG_A": 'R', 'DRUG_B': 'R'}

    assert test_catalogue.predict("M2@G74!:4&M2@G74X:2") == {"DRUG_B": "U"}
    assert test_catalogue2.predict("M2@G74!:0.04&M2@G74X:0.02") == {"DRUG_B": "U"}

def test_misc():
    '''Testing various edge cases
    '''
    test_catalogue = piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv")

    #Genes not in the catalogue should default to 'S'
    assert test_catalogue.predict("K@g-7a") == 'S'
    assert test_catalogue.predict("K@g-7a&K@S6L") == 'S'

    #Our catalogue doesn't have default rules for 'N'
    #So should cause errors
    with pytest.raises(ValueError):
        test_catalogue.predict("N@K45L")
    
    #Malformed minor mutation should cause issues
    with pytest.raises(AssertionError):
        test_catalogue.predict("M2@F75V:nope")

    #Trying with catalogues with issues, each hitting a different ValueError
    with pytest.raises(ValueError):
        piezo.ResistanceCatalogue("tests/test-catalogue/broken-catalogue.csv")

    with pytest.raises(ValueError):
        piezo.ResistanceCatalogue("tests/test-catalogue/broken-catalogue2.csv")

    with pytest.raises(ValueError):
        piezo.ResistanceCatalogue("tests/test-catalogue/broken-catalogue3.csv")

    with pytest.raises(ValueError):
        piezo.ResistanceCatalogue("tests/test-catalogue/broken-catalogue4.csv")
    
    with pytest.raises(ValueError):
        piezo.ResistanceCatalogue("tests/test-catalogue/broken-catalogue5.csv")

    with pytest.raises(ValueError):
        piezo.ResistanceCatalogue("tests/test-catalogue/broken-catalogue6.csv")
    
    #Checking for deprication warnings
    with pytest.warns(UserWarning):
        test_catalogue.predict("M2@L73P", verbose=True)
    
    #Minor edge cases of large dels
    test_catalogue = piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS-minor-COV.csv")
    assert test_catalogue.predict("M2@del_1.0") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.9:3") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.85:10") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.8:1") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.7") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.5") == {"DRUG_A": "U", "DRUG_B": "U"}

    test_catalogue = piezo.ResistanceCatalogue("tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS-minor-FRS.csv")
    assert test_catalogue.predict("M2@del_1.0") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.9:0.3") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.85:0.99") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.8:0.01") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.7") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.5:0.03") == {"DRUG_A": "U", "DRUG_B": "U"}


