#! /usr/bin/python3.5

import logging
import argparse

import piezo


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogue",default=None,required=False,help="the path to the resistance catalogue")
    parser.add_argument("--mutation",default=None,required=False,help="the mutation, in the correct form to match the grammar of the catalogue")
    parser.add_argument("--resistant_genes_only",action='store_true',default=False,help="only consider genes in the catalogue which have at least one variant associated with resistance")
    options = parser.parse_args()

    # instantiate a Resistance Catalogue instance by passing a text file
    resistance_catalogue=piezo.ResistanceCatalogue(options.catalogue,prediction_subset_only=options.resistant_genes_only)

    result_dict=resistance_catalogue.predict(options.mutation)

    print(result_dict)
