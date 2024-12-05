from collections import defaultdict
import pytest

import piezo


@pytest.mark.parametrize(
    "test_catalogue,genes,gene_lookup",
    [
        (
            piezo.ResistanceCatalogue(
                "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"
            ),
            ["M1", "M2", "N", "MULTI"],
            {
                "M1": ["DRUG_A", "DRUG_B"],
                "M2": ["DRUG_A", "DRUG_B"],
                "MULTI": ["DRUG_A", "DRUG_B"],
                "N": ["DRUG_A"],
            },
        ),
        (
            piezo.ResistanceCatalogue(
                "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv",
                prediction_subset_only=True,
            ),
            ["M1", "M2", "MULTI"],
            {
                "M1": ["DRUG_A", "DRUG_B"],
                "M2": ["DRUG_A", "DRUG_B"],
                "MULTI": ["DRUG_A", "DRUG_B"],
            },
        ),
    ],
)
def test_catalogue__init__(test_catalogue, genes, gene_lookup):
    assert test_catalogue.catalogue.genbank_reference == "NC_004148.2"

    assert test_catalogue.catalogue.name == "TEST"

    assert test_catalogue.catalogue.version == "v1.0"

    assert test_catalogue.catalogue.values == ["R", "F", "U", "S"]

    assert test_catalogue.catalogue.grammar == "GARC1"

    assert test_catalogue.catalogue.number_rows == 56

    assert test_catalogue.catalogue.drugs == ["DRUG_A", "DRUG_B"]

    # Ordering of lists is annoying here, so sort both
    genes = sorted(genes)
    assert sorted(test_catalogue.catalogue.genes) == genes

    assert test_catalogue.catalogue.gene_lookup == gene_lookup

    # Convert the gene_lookup to drug_lookup
    lookup = defaultdict(list)
    for gene in gene_lookup.keys():
        for drug in gene_lookup[gene]:
            lookup[drug].append(gene)
    lookup = {drug: sorted(lookup[drug]) for drug in lookup.keys()}

    actual_lookup = {
        drug: sorted(test_catalogue.catalogue.drug_lookup[drug])
        for drug in test_catalogue.catalogue.drug_lookup.keys()
    }
    assert actual_lookup == lookup


@pytest.mark.parametrize(
    "test_catalogue",
    [
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"
        ),
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv",
            prediction_subset_only=True,
        ),
    ],
)
def test_catalogue_prediction_snps(test_catalogue):
    # check a row in the catalogue
    assert test_catalogue.predict("M2@L73P") == {"DRUG_A": "U", "DRUG_B": "U"}

    # check a synonymous mutation has no effect
    assert test_catalogue.predict("M2@L73L") == {"DRUG_A": "S", "DRUG_B": "S"}

    # check a het
    assert test_catalogue.predict("M2@L73Z") == {"DRUG_A": "F", "DRUG_B": "S"}

    # check a filter fail that will hit the default S rule
    assert test_catalogue.predict("M2@L73O") == {"DRUG_A": "S", "DRUG_B": "S"}

    # check a null that will hit the default S rule
    assert test_catalogue.predict("M2@L73X") == {"DRUG_A": "S", "DRUG_B": "S"}

    # check a filter fail that hits a specific rule for DRUG_A
    assert test_catalogue.predict("M2@G74O") == {"DRUG_A": "F", "DRUG_B": "S"}

    # check a null that hits a specific rule for DRUG_A
    assert test_catalogue.predict("M2@G74X") == {"DRUG_A": "F", "DRUG_B": "S"}

    # If a null doesn't hit a specific rule, it shouldn't hit defaults of U
    assert test_catalogue.predict("M1@A12X") == "S"
    assert test_catalogue.predict("M1@A12X", show_evidence=True) == "S"

    # check hitting a wildtype row
    assert test_catalogue.predict("M2@G74I") == {"DRUG_A": "U", "DRUG_B": "U"}

    assert test_catalogue.predict("M2@G74Z") == {"DRUG_A": "R", "DRUG_B": "S"}

    assert test_catalogue.predict("M2@G74P") == {"DRUG_A": "R", "DRUG_B": "U"}

    assert test_catalogue.predict("M2@G74!") == {"DRUG_A": "S", "DRUG_B": "U"}

    assert test_catalogue.predict("M2@t-15c") == {"DRUG_A": "U", "DRUG_B": "U"}

    assert test_catalogue.predict("M2@t-15z") == {"DRUG_A": "S", "DRUG_B": "S"}

    assert test_catalogue.predict("M2@t-15o") == {"DRUG_A": "F", "DRUG_B": "S"}

    assert test_catalogue.predict("M2@t-15x") == {"DRUG_A": "F", "DRUG_B": "S"}

    # check that a gene not in the catalogue simply returns an "S"
    assert test_catalogue.predict("N1@S2T") == "S"

    # Checking large deletions
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
        assert test_catalogue.predict("M2@L73P") == {"DRUG_A": "R", "DRUG_B": "R"}

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


@pytest.mark.parametrize(
    "test_catalogue",
    [
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"
        ),
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv",
            prediction_subset_only=True,
        ),
    ],
)
def test_catalogue_prediction_snps_evidence(test_catalogue):
    # check a row in the catalogue
    assert test_catalogue.predict("M2@L73P", show_evidence=True) == {
        "DRUG_A": ("U", {"row": 12}),
        "DRUG_B": ("U", {"row": 25}),
    }

    # check a synonymous mutation has no effect
    assert test_catalogue.predict("M2@L73L", show_evidence=True) == {
        "DRUG_A": ("S", {"row": 11}),
        "DRUG_B": ("S", {"row": 24}),
    }

    # check a het
    assert test_catalogue.predict("M2@L73Z", show_evidence=True) == {
        "DRUG_A": ("F", {"row": 1}),
        "DRUG_B": ("S", {"row": 26}),
    }

    # check a filter fail that will hit the default S rule
    assert test_catalogue.predict("M2@L73O", show_evidence=True) == {
        "DRUG_A": ("S", {"row": 15}),
        "DRUG_B": ("S", {"row": 27}),
    }

    # check a null that will hit the default S rule
    assert test_catalogue.predict("M2@L73X", show_evidence=True) == {
        "DRUG_A": ("S", {"row": 16}),
        "DRUG_B": ("S", {"row": 28}),
    }

    # check a filter fail that hits a specific rule for DRUG_A
    assert test_catalogue.predict("M2@G74O", show_evidence=True) == {
        "DRUG_A": ("F", {"row": 5}),
        "DRUG_B": ("S", {"row": 27}),
    }

    # check a null that hits a specific rule for DRUG_A
    assert test_catalogue.predict("M2@G74X", show_evidence=True) == {
        "DRUG_A": ("F", {"row": 4}),
        "DRUG_B": ("S", {"row": 28}),
    }

    # check hitting a wildtype row
    assert test_catalogue.predict("M2@G74I", show_evidence=True) == {
        "DRUG_A": ("U", {"row": 6}),
        "DRUG_B": ("U", {"row": 25}),
    }

    assert test_catalogue.predict("M2@G74Z", show_evidence=True) == {
        "DRUG_A": ("R", {"row": 3}),
        "DRUG_B": ("S", {"row": 26}),
    }

    assert test_catalogue.predict("M2@G74P", show_evidence=True) == {
        "DRUG_A": ("R", {"row": 2}),
        "DRUG_B": ("U", {"row": 25}),
    }

    assert test_catalogue.predict("M2@G74!", show_evidence=True) == {
        "DRUG_A": ("S", {"row": 7}),
        "DRUG_B": ("U", {"row": 25}),
    }

    assert test_catalogue.predict("M2@t-15c", show_evidence=True) == {
        "DRUG_A": ("U", {"row": 17}),
        "DRUG_B": ("U", {"row": 29}),
    }

    assert test_catalogue.predict("M2@t-15z", show_evidence=True) == {
        "DRUG_A": ("S", {"row": 18}),
        "DRUG_B": ("S", {"row": 30}),
    }

    assert test_catalogue.predict("M2@t-15o", show_evidence=True) == {
        "DRUG_A": ("F", {"row": 10}),
        "DRUG_B": ("S", {"row": 31}),
    }

    assert test_catalogue.predict("M2@t-15x", show_evidence=True) == {
        "DRUG_A": ("F", {"row": 9}),
        "DRUG_B": ("S", {"row": 32}),
    }

    # check that a gene not in the catalogue simply returns an "S"
    assert test_catalogue.predict("N1@S2T", show_evidence=True) == "S"

    # Checking large deletions
    assert test_catalogue.predict("M2@del_1.0", show_evidence=True) == {
        "DRUG_A": ("R", {"row": 42}),
        "DRUG_B": ("U", {"row": 43}),
    }
    assert test_catalogue.predict("M2@del_0.9", show_evidence=True) == {
        "DRUG_A": ("R", {"row": 42}),
        "DRUG_B": ("U", {"row": 43}),
    }
    assert test_catalogue.predict("M2@del_0.85", show_evidence=True) == {
        "DRUG_A": ("R", {"row": 42}),
        "DRUG_B": ("U", {"row": 43}),
    }
    assert test_catalogue.predict("M2@del_0.8", show_evidence=True) == {
        "DRUG_A": ("R", {"row": 42}),
        "DRUG_B": ("U", {"row": 43}),
    }
    assert test_catalogue.predict("M2@del_0.7", show_evidence=True) == {
        "DRUG_A": ("U", {"row": 44}),
        "DRUG_B": ("U", {"row": 43}),
    }
    assert test_catalogue.predict("M2@del_0.5", show_evidence=True) == {
        "DRUG_A": ("U", {"row": 44}),
        "DRUG_B": ("U", {"row": 43}),
    }


@pytest.mark.parametrize(
    "test_catalogue",
    [
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"
        ),
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv",
            prediction_subset_only=True,
        ),
    ],
)
def test_catalogue_prediction_indels(test_catalogue):
    # check a general indel
    assert test_catalogue.predict("M2@300_indel") == {"DRUG_A": "U", "DRUG_B": "U"}

    # check an insertion
    assert test_catalogue.predict("M2@300_ins") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@-10_ins") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@-10_ins_2") == {"DRUG_A": "U", "DRUG_B": "U"}

    # check a deletion
    assert test_catalogue.predict("M2@300_del") == {"DRUG_A": "U", "DRUG_B": "U"}

    # check numbered insertions and deletions
    assert test_catalogue.predict("M2@300_ins_6") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@300_del_6") == {"DRUG_A": "U", "DRUG_B": "U"}

    # check a specific insertion
    assert test_catalogue.predict("M2@300_ins_acc") == {"DRUG_A": "U", "DRUG_B": "U"}

    assert test_catalogue.predict("M2@100_indel") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@100_ins") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@100_ins_3") == {"DRUG_A": "R", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@100_ins_act") == {"DRUG_A": "R", "DRUG_B": "U"}

    # check these hit the frameshift rule
    assert test_catalogue.predict("M2@100_ins_1") == {"DRUG_A": "R", "DRUG_B": "R"}
    assert test_catalogue.predict("M2@100_ins_acta") == {"DRUG_A": "R", "DRUG_B": "R"}
    assert test_catalogue.predict("M2@100_ins_400") == {"DRUG_A": "R", "DRUG_B": "R"}
    assert test_catalogue.predict("M2@100_ins_actgactg") == {
        "DRUG_A": "R",
        "DRUG_B": "R",
    }
    assert test_catalogue.predict("M2@300_ins_61") == {"DRUG_A": "U", "DRUG_B": "R"}
    assert test_catalogue.predict("M2@300_del_65") == {"DRUG_A": "U", "DRUG_B": "R"}
    assert test_catalogue.predict("M2@300_ins_acctt") == {"DRUG_A": "U", "DRUG_B": "R"}

    # badly formed INDELs
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@300_indel_5") == {
            "DRUG_A": "U",
            "DRUG_B": "U",
        }
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@300_indel_-5") == {
            "DRUG_A": "U",
            "DRUG_B": "U",
        }
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@300_insdel") == {"DRUG_A": "U", "DRUG_B": "U"}
    with pytest.raises(Exception):
        assert test_catalogue.predict("M2@300_del_ACT") == {
            "DRUG_A": "U",
            "DRUG_B": "U",
        }


@pytest.mark.parametrize(
    "test_catalogue",
    [
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"
        ),
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv",
            prediction_subset_only=True,
        ),
    ],
)
def test_multi(test_catalogue):
    # Exact match with 2 lines of the catalogue
    assert test_catalogue.predict("M2@G74!&M2@G74X") == {"DRUG_A": "U", "DRUG_B": "R"}

    # Should hit an 'S' (and 'U' default) and an 'R'
    assert test_catalogue.predict("M2@F75T&M2@G74P") == {"DRUG_A": "R", "DRUG_B": "U"}

    # Should hit an 'F' and an 'R' (and 'U' default)
    assert test_catalogue.predict("M2@L73Z&M2@G74P") == {"DRUG_A": "R", "DRUG_B": "U"}

    # Shouldn't hit anything other than 'S'
    assert test_catalogue.predict("M2@A5A&M2@K6K") == "S"

    # Should hit a general `R` rule for DRUG_A and a general epistasis rule for DRUG_B
    # There's also a specific multi for `R` for DRUG_B, so directly checks that epistasis > R
    assert test_catalogue.predict("M2@142_del&M2@M1L") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@142_del&M2@M1L", show_evidence=True) == {
        "DRUG_A": ("R", {"row": 45}),
        "DRUG_B": ("S", {"row": 47}),
    }


@pytest.mark.parametrize(
    "test_catalogue",
    [
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"
        ),
        piezo.ResistanceCatalogue(
            "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv",
            prediction_subset_only=True,
        ),
    ],
)
def test_minor_population(test_catalogue):
    # COV
    # Exact match
    assert test_catalogue.predict("M2@F75V:2") == {"DRUG_A": "R", "DRUG_B": "S"}

    # Greater than
    assert test_catalogue.predict("M2@F75V:3") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@F75V:4") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@F75V:350") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@F75V:123456789123456789") == {
        "DRUG_A": "R",
        "DRUG_B": "S",
    }

    # Less than (so should just hit a default for now)
    assert test_catalogue.predict("M2@F75V:1") == {"DRUG_A": "S", "DRUG_B": "S"}

    # Should fail
    with pytest.raises(Exception):
        test_catalogue.predict("M2@F75V:0")
    with pytest.raises(Exception):
        test_catalogue.predict("M2@F75V:-1")

    # Multi
    assert test_catalogue.predict("M2@G74!:4&M2@G74X:2") == {
        "DRUG_B": "U",
    }
    assert test_catalogue.predict("M2@G74!:17&M2@G74X:2") == {
        "DRUG_B": "U",
    }

    # FRS
    test_catalogue2 = piezo.ResistanceCatalogue(
        "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS-FRS.csv"
    )

    # Exact match
    assert test_catalogue2.predict("M2@F75V:0.03") == {"DRUG_A": "R", "DRUG_B": "S"}

    # Greater than
    assert test_catalogue2.predict("M2@F75V:0.04") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue2.predict("M2@F75V:0.05") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue2.predict("M2@F75V:0.9") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue2.predict("M2@F75V:0.123456789123456789") == {
        "DRUG_A": "R",
        "DRUG_B": "S",
    }

    # Less than (so should just hit a default for now)
    assert test_catalogue2.predict("M2@F75V:0.01") == {"DRUG_A": "S", "DRUG_B": "S"}
    assert test_catalogue2.predict("M2@F75V:0.02") == {"DRUG_A": "S", "DRUG_B": "S"}

    # Should fail
    with pytest.raises(Exception):
        test_catalogue2.predict("M2@F75V:0.0")
    with pytest.raises(Exception):
        test_catalogue2.predict("M2@F75V:-0.1")

    assert test_catalogue2.predict("M2@G74!:0.04&M2@G74X:0.02") == {
        "DRUG_B": "U",
    }
    assert test_catalogue2.predict("M2@G74!:0.84&M2@G74X:0.20") == {
        "DRUG_B": "U",
    }

    assert test_catalogue.predict("M2@G74!:4&M2@G74X:2") == {"DRUG_B": "U"}
    assert test_catalogue2.predict("M2@G74!:0.04&M2@G74X:0.02") == {"DRUG_B": "U"}


def test_misc():
    """Testing various edge cases"""
    test_catalogue = piezo.ResistanceCatalogue(
        "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv"
    )

    # Genes not in the catalogue should default to 'S'
    assert test_catalogue.predict("K@g-7a") == "S"
    assert test_catalogue.predict("K@g-7a&K@S6L") == "S"

    # Our catalogue doesn't have default rules for 'N'
    # So should cause errors
    with pytest.raises(ValueError):
        test_catalogue.predict("N@K45L")

    # Malformed minor mutation should cause issues
    with pytest.raises(AssertionError):
        test_catalogue.predict("M2@F75V:nope")

    # Trying with catalogues with issues, each hitting a different ValueError
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

    with pytest.raises(ValueError) as e:
        piezo.ResistanceCatalogue("tests/test-catalogue/broken-catalogue7.csv")
    assert (
        str(e.value)
        == "Badly formed mutation: M2@*?&M2@S450L contains generic rules which cover specific rules!"
    )

    # Checking for deprication warnings
    with pytest.warns(UserWarning):
        test_catalogue.predict("M2@L73P", verbose=True)

    # Minor edge cases of large dels
    test_catalogue = piezo.ResistanceCatalogue(
        "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS-minor-COV.csv"
    )
    assert test_catalogue.predict("M2@del_1.0") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.9:3") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@del_0.85:10") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@del_0.8:1") == {"DRUG_A": "S", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@del_0.7") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.5") == {"DRUG_A": "U", "DRUG_B": "U"}

    test_catalogue = piezo.ResistanceCatalogue(
        "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS-minor-FRS.csv"
    )
    assert test_catalogue.predict("M2@del_1.0") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.9:0.3") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@del_0.85:0.99") == {"DRUG_A": "R", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@del_0.8:0.01") == {"DRUG_A": "S", "DRUG_B": "S"}
    assert test_catalogue.predict("M2@del_0.7") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@del_0.5:0.03") == {"DRUG_A": "S", "DRUG_B": "S"}

    # Ensure that if a deletion crosses the boundary into CDS, it should predict
    # the coding part too. i.e if a deletion starts in the promoter, but causes a frameshift
    # in the coding region, it should predict the frameshift

    # This catalogue has a frameshift rule for M2, so this should predict 'R' for DRUG_B
    # in cases the deletion crosses the boundary

    # Just promoter
    assert test_catalogue.predict("M2@-12_del_2") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@-12_del_7") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@-5_del_acgt") == {"DRUG_A": "U", "DRUG_B": "U"}

    # Deleting all of promoter (shouldn't hit frameshift)
    assert test_catalogue.predict("M2@-12_del_12") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@-1_del_a") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@-10_del_" + "a" * 10) == {
        "DRUG_A": "U",
        "DRUG_B": "U",
    }

    # Deletion passing into CDS
    assert test_catalogue.predict("M2@-12_del_13") == {"DRUG_A": "U", "DRUG_B": "R"}
    assert test_catalogue.predict("M2@-12_del_14") == {"DRUG_A": "U", "DRUG_B": "R"}
    assert test_catalogue.predict("M2@-5_del_15") == {"DRUG_A": "U", "DRUG_B": "R"}
    assert test_catalogue.predict("M2@-5_del_acgtac") == {"DRUG_A": "U", "DRUG_B": "R"}
    # Not a frameshift, so should remain 'U'
    assert test_catalogue.predict("M2@-12_del_15") == {"DRUG_A": "U", "DRUG_B": "U"}
    assert test_catalogue.predict("M2@-3_del_12") == {"DRUG_A": "U", "DRUG_B": "U"}

    # Specific case of a deletion that hits a specific rule for DRUG_A as well as a frameshift
    assert test_catalogue.predict("M2@-5_del_10") == {"DRUG_A": "R", "DRUG_B": "R"}

    # Wildcard ins shouldn't give prediction on del
    # (yes this sounds ridiculous but was present for years)
    # This catalogue has a dummy drug ("M3") which doesn't have default rules,
    #   so this should crash if it misses the `*_ins`
    test_catalogue = piezo.ResistanceCatalogue(
        "tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS-FRS.csv"
    )
    assert test_catalogue.predict("M3@12_ins_c") == {"DRUG_B": "R"}
    with pytest.raises(ValueError):
        print(test_catalogue.predict("M3@12_del_c"))

    # Double checking that a minor allele doesn't hit a general rule anymore
    assert test_catalogue.predict("M2@37_del_c:0.2") == {"DRUG_A": "S", "DRUG_B": "S"}
