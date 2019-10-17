[![Build Status](https://travis-ci.com/philipwfowler/piezo.svg?branch=develop)](https://travis-ci.com/philipwfowler/piezo) [![codecov](https://codecov.io/gh/philipwfowler/piezo/branch/develop/graph/badge.svg)](https://codecov.io/gh/philipwfowler/piezo) [![Documentation Status](https://readthedocs.org/projects/piezo/badge/?version=latest)](https://piezo.readthedocs.io/en/latest/?badge=latest)

# piezo

Takes a Clockwork VCF, a Genetic Catalogue and the associated GenBank file and predicts the effect on the relevant antibiotics. This code was developed as part of the [CRyPTIC](http://www.crypticproject.org) international tuberculosis consortium. If you would like to use the software commercially, please consult the LICENCE file.

## Example

A demonstration can be found in the piezo-vcf-parse.py script in the bin/ folder of the repository which following installation should be in your $PATH.

```
# setup a catalogue
cat=piezo.ResistanceCatalogue(input_file="config/NEJM2018-RSU-catalogue-H37rV_v3.csv",
                              genome_object="config/H37rV_v3.gbk",
                              catalogue_name="NEJM2018")

# use the predict method
print(cat.predict(gene_mutation='rpoB_S450L'))

{'RIF': 'R'}
```

## Pre-requisites

Everything is Python3. The only non-standard pre-requisite is [gumpy](https://github.com/philipwfowler/gumpy), which is on pypi so should install automagickally anywhere, but you could always clone and install it first... Previously the package relied on [snpit](https://github.com/philipwfowler/snpit), [datreant](https://datreant.readthedocs.io/en/latest/), [gemucator](https://github.com/philipwfowler/gemucator) and [vasta](https://github.com/philipwfowler/vasta). SNPIT is specific to TB so has been removed. Users may not wish to use Datreant (although I recommend it) so its functionality has gone into upstream code and the functionality of the last two has been rolled into gumpy.

Everything else is fairly standard / should be installed by the `setup.py` if you don't already have it
* pandas
* numpy
* h5py
* pysam
* biopython
* tqdm
* pytest
* pytest-cov


## Installation

First clone/download the repo. Then

```
$ cd piezo
$ ls
$ python setup.py develop --user
```

If you've done prerequisites 1 and 2, it hopefully will just work...

## Included files

```
$ ls config/
NC_004148.2.gbk                         LID2015-RSU-catalogue-H37rV_v2.csv
H37rV_v3.gbk                            NEJM2018-RSU-catalogue-H37rV_v3.csv
H37rV_v3.pkl.gz                         CRYPTICv1.0-RSU-catalogue-H37rV_v3.csv
NC_004148.2_RSU_catalogue.csv
```
NC_004148 is the reference genome of the human metapneumovirus and is used primarily for unit testing since it is small and fast to parse. The other catalogues are based on these papers

Walker TM, Kohl TA, Omar S V, Hedge J, Del Ojo Elias C, et al. Whole-genome sequencing for prediction of Mycobacterium tuberculosis drug susceptibility and resistance: a retrospective cohort study. Lancet Infec Dis 2015;15:1193–202. [https://doi.org/10.1016/S1473-3099(15)00062-6](doi:10.1016/S1473-3099(15)00062-6

The CRyPTIC Consortium, 100000 Genomes Project. Prediction of Susceptibility to First-Line Tuberculosis Drugs by DNA Sequencing. N Engl J Med 2018;379:1403–1415. [http://doi.org/10.1056/NEJMoa1800474](doi:10.1056/NEJMoa1800474)

are included in `config/`. The first is only relevant when used alongside version 2 of the H37rV GenBank catalogue and the second with version 3 (it will fail if you try using the wrong version since some of the genes are renamed and shifted). The CRYPTICv1.0 catalogue is an amalgam of the above two, with respect to version 3 of the H37rV TB reference genome.

2. H37rV GenBank files

Only version 3 is included which is the current version and was only released in Dec 2017. CRyPTIC is using version 3.

## End-to-end example

In the `piezo` folder there is are several VCF files in e.g. `examples/01/01.vcf` and `examples/02/02.vcf`. To analyse it run this command

```
$ piezo-vcf-parse.py  --vcf_file examples/01/01.vcf\
                      --genome_object config/H37rV_v3.pkl.gz\
                      --catalogue_file config/NEJM2018-RSU-catalogue-H37rV_v3.csv\
                      --catalogue_name NEJM2018\
                      --ignore_vcf_filter\
                      --ignore_vcf_status\
                      --progress
```

These two VCF files were produced with an early version of Clockwork which didn't use the `FILTER` and `PASS` fields in the VCF format, hence the flags.

You should see a series of lines written to STDOUT with the high-level predictions (In this case the sample is resistant to all four drugs in the NEJM2018 catalogue). In addition, two CSVs are written here

```
$ ls examples/01/
01-effects.csv   01-mutations.csv  01-variants.csv  01.vcf
```

The variants CSV has one row per single nucleotide polymorphism or insertion/deletion. The mutations CSV is one row per amino acid mutation, where appropriate, and therefore usually has fewer rows than the variants CSV, whereas the effects CSV is one row per mutation per drug with the associated predicted phenotype listed. Hence the later is usually a little longer than the former, since some mutations confer effects on multiple drugs.
