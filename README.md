[![Tests](https://github.com/oxfordmmm/piezo/actions/workflows/tests.yaml/badge.svg)](https://github.com/oxfordmmm/piezo/actions/workflows/tests.yaml)
[![codecov](https://codecov.io/gh/oxfordmmm/piezo/branch/master/graph/badge.svg)](https://codecov.io/gh/oxfordmmm/piezo)
[![Documentation Status](https://readthedocs.org/projects/piezo/badge/?version=latest)](https://piezo.readthedocs.io/en/latest/?badge=latest) [![PyPI version](https://badge.fury.io/py/piezo.svg)](https://badge.fury.io/py/piezo)

# piezo

Predict the effect of a genetic mutation on the effect of an antibiotic using a supplied AMR catalogue.

This code was developed as part of the [CRyPTIC](http://www.crypticproject.org) international tuberculosis consortium. If you would like to use the software commercially, please consult the LICENCE file.

## Installation

### using `pip`

This will install the most recent release on PyPI.

```
$ pip install piezo
$ cd piezo
$ py.test
```

### from GitHub

This will install the current version from GitHub and therefore may be ahead of the PyPI version.

```
$ git clone https://github.com/oxfordmmm/piezo
$ cd piezo
$ pip install .
$ py.test
```
The pre-requisites are all fairly standard and are listed in `setup.cfg` so will be automatically installed.

## Included files

```
$ $ ls tests/test-catalogue/
NC_004148.2.gbk                    NC_004148.2_TEST_GM1_RFUS_v1.0.csv
```
NC_004148 is the reference genome of the human metapneumovirus and is used primarily for unit testing since it is small and fast to parse.

## Design of AMR catalogue

`piezo` is written so as to be extendable in the future to other ways of describing genetic variation with respect to a reference. It includes the concept of a `grammar` which specifies how the genetic variation is described.

At present only a single grammar, `GARC1` is supported. `GARC` is short for Grammar for Antimicrobial Resistance Catalogues. This grammar is described in more detail [elsewhere](http://fowlerlab.org/2018/11/25/goarc-a-general-ontology-for-antimicrobial-resistance-catalogues/), however in brief, it is a gene-centric view (and therefore has no way of describing genetic variation that lies outside a coding region, other than as a 'promoter' mutation). All mutations start with the gene (or locus) name which must match the name of a gene (or locus) in the relevant GenBank file. It is the user's responsibility to ensure this, although e.g. the `gumpy` package can be used to perform such sanity checks. The mutation is delineated from the gene using a `@` symbol and within the mutation `_` is used as a field separator to separate the different components. All variation is described as either a `SNP` or an `INDEL`. If they occur within a coding region `SNP`s are specified by their effect on the amino acids which are always in UPPERCASE e.g. `rpoB@S450L`. If in the assumed promoter region, then the nucleotide change and position is specified e.g. `fabG1@c-15t`. Nucleotides are always in lowercase. `INDEL`s can be specified at different levels of granularity e.g. `rpoB@1250_indel` means 'any insertion of deletion at this position', but we could equally be highly specific and say `rpoB@1250_ins_cta` which means 'an insertion of cta at this position'. There is also the special case of frameshifting mutations which are described by `fs`.

Wildcards are also supported. Hence `rpoB@*?` means 'any non-synoymous mutation in the coding region of the protein'. To avoid confusion the stop codon is represented by `!` which is non-standard. Het calls are, at present, represented by a `Z` or `z` depending on whether they occur in the coding or promoter regions. This may be extended in the future. Likewise null calls are represented by an `X` or `x`.

The general principle is each mutation can 'hit' multiple rules in the catalogue, but it is the most specific rule that will be followed. Hence consider a toy example, again from TB

```
rpoB@*?     RIF   U   any non-synoymous mutation in the coding region has an unknown effect of RIF
rpoB@S450?  RIF   R   any non-synoymous mutation at Ser450 confers resistance
rpoB@S450Z  RIF   F   a het call at Ser450 should be reported as an F (fail).
```

## Example

A demonstration script called `piezo-predict.py` can be found in the `bin/` folder of the repository which following installation should be in your `$PATH`. A made-up catalogue for testing purposes can be found in `tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv` which is based on the Human metapneumovirus, however the entries are fictious. It contains two drugs and a series of mutations in the *M2* gene.

```
$ piezo-predict.py --catalogue tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv --mutation M2@L73L
{'DRUG_B': 'S', 'DRUG_A': 'S'}

$ piezo-predict.py --catalogue tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv --mutation M2@L73R
{'DRUG_A': 'R', 'DRUG_B': 'U'}

$ piezo-predict.py --catalogue tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv --mutation M2@L73Z
{'DRUG_B': 'S', 'DRUG_A': 'F'}

$ piezo-predict.py --catalogue tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv --mutation M2@300_indel
{'DRUG_B': 'U', 'DRUG_A': 'U'}

$ piezo-predict.py --catalogue tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv --mutation M2@300_ins
{'DRUG_B': 'U', 'DRUG_A': 'U'}

$ piezo-predict.py --catalogue tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv --mutation M2@300_ins_2
{'DRUG_B': 'U', 'DRUG_A': 'U'}

$ piezo-predict.py --catalogue tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv --mutation M2@300_ins_3
{'DRUG_A': 'U', 'DRUG_B': 'R'}

$ piezo-predict.py --catalogue tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv --mutation M2@300_ins_4
{'DRUG_B': 'U', 'DRUG_A': 'U'}

$ piezo-predict.py --catalogue tests/test-catalogue/NC_004148.2_TEST_v1.0_GARC1_RFUS.csv --mutation M2@300_ins_cta
{'DRUG_B': 'R', 'DRUG_A': 'U'}
```
