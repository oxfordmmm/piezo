# piezo

Takes a Clockwork VCF, a Genetic Catalogue and the associated GenBank file and predicts the effect on the relevant antibiotics.

## Example

A more detailed description can be found in the Jupyter Notebook in the root of the repository.

```
# setup a catalogue
cat=piezo.ResistanceCatalogue(input_file="config/LID2015-RSU-catalogue-v1.0-H37rV_v2.csv",
                              genbank_file="config/H37rV_v2.gbk",
                              catalogue_name="LID2015B")

# use the predict method
print(cat.predict(gene_mutation='rpoB_S450L'))

{'RIF': 'R'}
```

Then an end-to-end script is supplied in `bin/piezo-vcf-parse.py` which following installation should be in your $PATH.

## Pre-requisites

Everything is Python3.

1. [gemucator](https://github.com/philipwfowler/gemucator). This is not on pypi or anything, so you need to clone and install it first...

2. [snpit](https://github.com/philipwfowler/snpit). This is Class based version of Sam Lipworth's SNP-IT code. It is only used in the `bin/piezo-vcf-parse.py` example and can be easily removed. Again, you'll need to clone and install first.

3. [datreant](https://datreant.readthedocs.io/en/latest/). Is used heavily in CRyPTIC and not used so much here. I probably could cut it out, but I like it... Is on Pypi so should automatically install.

4. Everything else is fairly standard / should be installed by the `setup.py`
* pandas
* numpy
* PyVCF
* Biopython
* tqdm

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
H37rV_v2.gbk                            LID2015-RSU-catalogue-v1.0-H37rV_v2.csv
H37rV_v3.gbk
```

1. A catalogue based on this paper

Walker TM, Kohl TA, Omar S V, Hedge J, Del Ojo Elias C, et al. Whole-genome sequencing for prediction of Mycobacterium tuberculosis drug susceptibility and resistance: a retrospective cohort study. Lancet Infec Dis 2015;15:1193–202. doi:10.1016/S1473-3099(15)00062-6

is included in `config/LID2015-RSU-catalogue-v1.0-H37rV_v2.csv`. As suggested by the filename, it is only relevant when used alongside version 2 of the H37rV GenBank catalogue (it will fail if you try using v3 since some of the genes are renamed and shifted).

2. H37rV GenBank files

Both versions 2 and 3 are included. Version 3 is the current version, but was only released in Dec 2017. CRyPTIC is using version 3.

## End-to-end example

In the `piezo` folder there is a single VCF file in `examples/01/01.vcf`. To analyse it run this command

```
$ piezo-vcf-parse.py  --vcf_file examples/01/01.vcf 
                      --genbank_file config/H37rV_v2.gbk 
                      --resistance_catalogue config/LID2015-RSU-catalogue-v1.0-H37rV_v2.csv 
                      --catalogue_name LID2015B
```

You should see a series of lines written to STDOUT with the high-level predictions, some information about Lineage (gleaned using [snpit](https://github.com/philipwfowler/snpit)). In addition, two CSVs are written here

```
$ ls examples/01/
01-effects.csv   01-mutations.csv 01.vcf
```

The mutations CSV is one row per mutation found in the genes listed in the catalogue, whereas the effects CSV is one row per mutation per drug with the associated predicted phenotype listed. Hence the later is usually a little longer than the former, since some mutations confer effects on multiple drugs.




