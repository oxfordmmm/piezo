#! /usr/bin/env python
'''Instanciate the catalogue
'''

import collections
import os
import warnings

import pandas
import ujson

from piezo.grammar_GARC1 import predict_GARC1, process_catalogue_GARC1

# define the named tuple that will specify a resistance catalogue
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
                                      'number_rows',
                                      'rules' ]
                                    )


class ResistanceCatalogue:
    '''Resistance catalogue loading
    '''

    def __init__(self, catalogue_file: str, prediction_subset_only: bool=False):
        '''Construct a resistance catalogue

        Args:
            catalogue_file (str): Path to the catalogue file
            prediction_subset_only (bool, optional): Whether to use a subset of genes to only resistance genes. Defaults to False.
        '''
        self.catalogue = load_catalogue(catalogue_file, prediction_subset_only)

    def predict(self, mutation: str, verbose: bool=False) -> {str: str}:
        '''Make a prediction of a mutation's effects based on the catalogue

        Args:
            mutation (str): Mutation in GARC
            verbose (bool, optional): Whether to be verbose. Defaults to False. DEPRECIATED

        Returns:
            {str: str}: Dictionary mapping drug name -> prediction, or "S" if no matches
        '''
        if verbose:
            warnings.warn("`verbose` kwarg is depreciated and will be removed in future.", UserWarning)
        return predict(self.catalogue, mutation=mutation, verbose=verbose)

def parse_json(data: str) -> dict:
    '''Load the data within a json string to Python dict

    Args:
        data (str): JSON string

    Returns:
        dict: Dictionary of the JSON data
    '''

    return ujson.loads(data)


def load_catalogue(catalogue_file: str, prediction_subset_only: bool) -> collections.namedtuple:
    '''
    Read in the Antimicrobial Resistance Catalogue.

    Args:
        catalogue_file (str): path to a resistance catalogue as a CSV file in the correct format, as determined by its grammar
        prediction_subset_only (bool): whether to subset the catalogue down so it ONLY includes entities (e.g. DRUG,GENE pairs) that include at least one row predicting resistance

    Returns:
        catalogue (collections.namedtuple): defined tuple

    Notes:
        * Applies checks to ensure the catalogue is in the right format. TODO: Remove assert statements as they should not be used for runtime error checking
    '''


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

    number_rows=len(rules)

    assert len(rules.PREDICTION_VALUES.unique())==1, "multiple grammar value types used in the catalogue!"
    values=rules.PREDICTION_VALUES.unique()[0]
    if values in ['RS','RUS','RFUS']:
        values=[i for i in rules.PREDICTION_VALUES.unique()[0]]
        assert sorted(values)==sorted(list(rules.PREDICTION.unique())), "PREDICTION column contains entries not in "+rules.PREDICTION_VALUES.unique()[0]
    else:
        raise ValueError("content of column CATALOGUE_VALUES not recognised!")


    # since we are storing the above values in the named tuple, we don't need them in the rules dataframe
    rules.drop(columns=['GENBANK_REFERENCE', 'CATALOGUE_NAME', 'CATALOGUE_VERSION', 'CATALOGUE_GRAMMAR', 'PREDICTION_VALUES'],inplace=True)

    # find out the DRUGS that this catalogue applies to
    drugs=list(rules.DRUG.unique())

    if grammar!="GARC1":

        raise ValueError("only the GENE_MUTATION grammar is supported at present!")

    elif grammar=="GARC1":

        (rules,genes,drug_lookup,gene_lookup)=process_catalogue_GARC1(rules,drugs,catalogue_genes_only=prediction_subset_only)

    return catalogue(genbank_reference,
                           name,
                           version,
                           grammar,
                           values,
                           drugs,
                           genes,
                           drug_lookup,
                           gene_lookup,
                           number_rows,
                           rules)

def predict(catalogue: collections.namedtuple, mutation: str, verbose: bool=False) -> {str: str}:
    '''
    Predict the effect of the given mutation on one or more antimicrobials.

    Args:
        catalogue (collections.namedtuple): The catalogue
        mutation (str): a genetic variant in the form GENE_MUTATION e.g. for a SNP katG_S315T, INDEL katG_315_indel.
        verbose (bool): if True, then a description of the rules that apply to the supplied mutation and their priority is written to STDOUT (default=False). DEPRECIATED

    Returns:
        result (dict): the drugs affected by the mutation are the keys, and the predicted phenotypes are the values. e.g. {'LEV':'R', 'MXF':'R'}
                       if the gene isn't in the catalogue, then an "S" is returned, on the assumption that it is probably susceptible.

    Notes:
        * mutations can be specified in a grammar that covers most of the known and expected genetic variants.
        * Stop codon is represented by "!"
        * "any mutation at position S315" (i.e. a wildcard) is represented by "?" e.g. S315?
        * for more info see the walkthrough and also the NOMENCLATURE.md file
    '''
    if verbose:
        warnings.warn("`verbose` kwarg is depreciated and will be removed in future.", UserWarning)
    
    return predict_GARC1(catalogue, mutation)
