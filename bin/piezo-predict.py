#! /usr/bin/python3.5

import logging
import argparse

import piezo


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--catalogue",default=None,required=False,help="the path to the resistance catalogue")
    parser.add_argument("--mutation",default=None,required=False,help="the mutation, in the correct form to match the grammar of the catalogue")
    options = parser.parse_args()

    # instantiate a Resistance Catalogue instance by passing a text file
    resistance_catalogue=piezo.ResistanceCatalogue(options.catalogue)

    result_dict=resistance_catalogue.predict(options.mutation)

    print(result_dict)
