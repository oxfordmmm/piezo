#! /usr/bin/python3.5

import argparse
from copy import deepcopy

import pandas
from tqdm import tqdm

import cryptic.genetics

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--database",default='CRyPTIC2',help="the name of the database (default='CRyPTIC2')")
    parser.add_argument("--study",default='V2',help="the name of the study (default='V2')")
    parser.add_argument("--input",required=True,help="either the path to a single VCF file, or a text file containing a list of paths to multiple VCF files to read in")
    parser.add_argument("--genbank",default="H37Rv.gbk",help="the genbank file of the H37Rv M. tuberculosis reference_collection genome")
    parser.add_argument("--resistance_catalogue",default="Catalogue_of_variants_3.csv",required=False,help="the path to the resistance catalogue")
    parser.add_argument("--verbose",action='store_true',default=False,help="the path of the output CSV file containing the list of all the detected mutations")
    options = parser.parse_args()

    assert options.database in ['CRyPTIC1','CRyPTIC2'], "invalid database!"
    assert options.study in ['V1','V2'], "invalid study"

    # setup a GeneCollection object that contains all the genes/loci we are interested in
    reference_collection=cryptic.genetics.GeneCollection(species="M. tuberculosis",genbank_file=options.genbank,database=options.database,study=options.study)

    # instantiate a Resistance Catalogue instance by passing a text file
    walker_catalogue=cryptic.genetics.ResistanceCatalogue(input_file=options.resistance_catalogue)

    list_of_vcf_files=[]

    # test if the filename we have been given ends in .vcf
    if options.input.split('.')[-1]=="vcf":

        # if it does, that means we are only going to process a single VCF file
        list_of_vcf_files.append(options.input)

    # otherwise assume we've been given a text file where each line is a relative path to a VCF file
    else:

        INPUT=open(options.input,'r')
        for line in INPUT:
            list_of_vcf_files.append(line.rstrip())
        INPUT.close()

    # count the number of vcf files in the specifed file
    number_vcf_files=len(list_of_vcf_files)

    # now that we have setup a GeneCollection of reference genes, we can work through the VCF files one-by-one
    for vcf_file in tqdm(list_of_vcf_files,total=number_vcf_files,disable=options.verbose):

        MUTATIONS_dict={}
        MUTATIONS_counter=0
        EFFECTS_dict={}
        EFFECTS_counter=0

        # create a copy of the reference_collection genes which we will then alter according to the VCF file
        # need a deepcopy to ensure we take all the private variables etc with us, and point just take pointers
        vcf_collection=deepcopy(reference_collection)

        vcf_collection.apply_vcf_file(vcf_file)

        (study,instance,UniqueID,site_id)=(vcf_collection.study,vcf_collection.instance,vcf_collection.UniqueID,vcf_collection.site_id)

        # by default assume wildtype behaviour so set all drug phenotypes to be susceptible
        phenotype={}
        for drug in walker_catalogue.drug_list:
            phenotype[drug]="S"

        # now get all the genes to calculate their own differences w.r.t the references, i.e. their mutations!
        for gene_name in reference_collection.gene_panel:

            # now we pass the reference_collection gene so our VCF gene can calculate its mutations
            vcf_collection.gene[gene_name].identify_mutations(reference_collection.gene[gene_name])

            mutations=vcf_collection.gene[gene_name].mutations

            for mutation_name in mutations:

                # store
                MUTATIONS_dict[MUTATIONS_counter]=[ UniqueID,site_id,gene_name,\
                                                    mutations[mutation_name]["TYPE"],\
                                                    mutation_name,\
                                                    mutations[mutation_name]["ELEMENT_TYPE"],\
                                                    mutations[mutation_name]["POSITION"],\
                                                    mutations[mutation_name]["PROMOTER"],\
                                                    mutations[mutation_name]["CDS"],\
                                                    mutations[mutation_name]["SYNONYMOUS"],\
                                                    mutations[mutation_name]["NONSYNONYMOUS"],\
                                                    mutations[mutation_name]["INSERTION"],\
                                                    mutations[mutation_name]["DELETION"],\
                                                    mutations[mutation_name]["REF"],\
                                                    mutations[mutation_name]["ALT"],\
                                                    mutations[mutation_name]["NUMBER_NUCLEOTIDE_CHANGES"]]
                MUTATIONS_counter+=1


                prediction=walker_catalogue.predict(mutation=gene_name+"_"+mutation_name)

                # if it isn't an S, then a dictionary must have been returned
                if prediction!="S":

                    # iterate through the drugs in the dictionary (can be just one)
                    for drug_name in prediction:

                        # only for completeness as this logic never leads to a change since by default the phenotype is S
                        if drug_name and prediction[drug_name]=="S" and phenotype[drug_name]=="S":
                            phenotype[drug_name]="S"

                        # if the prediction is a U, we only move to a U if the current prediction is S
                        # (to stop it overiding an R)
                        # (again for completeness including the superfluous state)
                        elif drug_name and prediction[drug_name]=="U" and phenotype[drug_name] in ["S","U"]:
                            phenotype[drug_name]="U"

                        # finally if an R is predicted, it must be R
                        elif drug_name and prediction[drug_name]=="R":
                            phenotype[drug_name]="R"

                        EFFECTS_dict[EFFECTS_counter]=[UniqueID,site_id,gene_name,mutation_name,drug_name,prediction[drug_name]]
                        EFFECTS_counter+=1
                else:
                    EFFECTS_dict[EFFECTS_counter]=[UniqueID,site_id,gene_name,mutation_name,"UNK","U"]
                    EFFECTS_counter+=1

        MUTATIONS=pandas.DataFrame.from_dict(MUTATIONS_dict,orient="index",columns=["UNIQUEID","SITEID","GENE","MUTATION_TYPE","MUTATION","ELEMENT_TYPE","POSITION","PROMOTER","CDS","SYNONYMOUS","NONSYNONYMOUS","INSERTION","DELETION","REF","ALT","NUMBER_NUCLEOTIDE_CHANGES"])
        MUTATIONS.set_index(["UNIQUEID","GENE","MUTATION"],inplace=True)
        MUTATIONS.to_csv(vcf_collection.vcf_folder+"/mutations.csv")

        EFFECTS=pandas.DataFrame.from_dict(EFFECTS_dict,orient="index",columns=["UNIQUEID","SITEID","GENE","MUTATION","DRUG","PREDICTION"])
        EFFECTS.set_index(["UNIQUEID","DRUG","GENE","MUTATION"],inplace=True)
        EFFECTS.to_csv(vcf_collection.vcf_folder+"/effects.csv")

        wgs_prediction_string=""
        for drug in walker_catalogue.drug_list:
            vcf_collection.datreant_file.categories["WGS_PREDICTION_"+drug]=phenotype[drug]
            wgs_prediction_string+=phenotype[drug]
        vcf_collection.datreant_file.categories["WGS_PREDICTION_STRING"]=wgs_prediction_string

        vcf_collection.determine_mdr()
