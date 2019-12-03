#! /usr/bin/env python

import os
import collections

import pandas, numpy, ujson

from piezo.gm1 import predict_GM1
from piezo.gm1 import process_catalogue_GM1


catalogue = collections.namedtuple('catalogue',
                                    [ 'genbank_reference',
                                      'name',
                                      'version',
                                      'grammar',
                                      'values',
                                      'drugs',
                                      'genes',
                                      'drug_lookup',
                                      'gene_lookup',
                                      'rules' ]
                                    )

class ResistanceCatalogue:

    def __init__(self, catalogue_file):

        self.catalogue = load_catalogue(catalogue_file)

    def predict(self, mutation, verbose=False):

        return predict(self.catalogue, mutation=mutation, verbose=verbose)

def parse_json(data):

    return ujson.loads(data)


def load_catalogue(catalogue_file):

    assert os.path.isfile(catalogue_file), "supplied catalogue file "+catalogue_file+" does not exist!"

    rules=pandas.read_csv(catalogue_file,converters={'OTHER':parse_json,'EVIDENCE':parse_json,'SOURCE':parse_json})

    assert len(rules.GENBANK_REFERENCE.unique())==1, "multiple genbank references specified in catalogue!"
    genbank_reference=rules.GENBANK_REFERENCE.unique()[0]

    assert len(rules.CATALOGUE_NAME.unique())==1, "multiple catalogue models found in the catalogue!"
    name=rules.CATALOGUE_NAME.unique()[0]

    assert len(rules.CATALOGUE_VERSION.unique())==1, "multiple catalogue versions found in the catalogue!"
    version=rules.CATALOGUE_VERSION.unique()[0]

    assert len(rules.CATALOGUE_GRAMMAR.unique())==1, "multiple grammars found in the catalogue!"
    grammar=rules.CATALOGUE_GRAMMAR.unique()[0]

    assert len(rules.CATALOGUE_VALUES.unique())==1, "multiple grammar value types used in the catalogue!"
    values=rules.CATALOGUE_VALUES.unique()[0]
    if values in ['RS','RUS','RFUS']:
        values=[i for i in rules.CATALOGUE_VALUES.unique()[0]]
    else:
        raise ValueError("content of column CATALOGUE_VALUES not recognised!")

    # since we are storing the above values in the named tuple, we don't need them in the rules dataframe
    rules=rules[['DRUG','MUTATION','PREDICTION','EVIDENCE','SOURCE','OTHER']]

    # find out the DRUGS that this catalogue applies to
    drugs=list(rules.DRUG.unique())

    if grammar!="GM1":

        raise ValueError("only the GENE_MUTATION grammar is supported at present!")

    elif grammar=="GM1":

        (rules,genes,drug_lookup,gene_lookup)=process_catalogue_GM1(rules,drugs)

        return catalogue(genbank_reference,
                               name,
                               version,
                               grammar,
                               values,
                               drugs,
                               genes,
                               drug_lookup,
                               gene_lookup,
                               rules)

def predict(catalogue,mutation,verbose):

    if catalogue.grammar=="GM1":

        result = predict_GM1(catalogue,mutation,verbose)

        return(result)

    else:

        raise ValueError("Only the GENE_MUTATION grammar is supported at present")
