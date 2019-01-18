#! /usr/bin/python3.5

import argparse,os
from copy import deepcopy

import pandas
from tqdm import tqdm

import piezo

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--database",default='CRyPTIC2',help="the name of the database (default='CRyPTIC2')")
    parser.add_argument("--study",default='V2',help="the name of the study (default='V2')")
    parser.add_argument("--vcf_file",required=True,help="the path to a single VCF file")
    parser.add_argument("--genbank",default="H37Rv.gbk",help="the genbank file of the H37Rv M. tuberculosis wildtype_gene_collection genome")
    parser.add_argument("--resistance_catalogue",default="Catalogue_of_variants_3_mod.csv",required=False,help="the path to the resistance catalogue")
    parser.add_argument("--verbose",action='store_true',default=False,help="whether to show progress using tqdm")
    options = parser.parse_args()

    # instantiate a Resistance Catalogue instance by passing a text file
    walker_catalogue=piezo.ResistanceCatalogue(input_file=options.resistance_catalogue)

    # retrieve the dictionary of genes from the Resistance Catalogue
    gene_panel=walker_catalogue.gene_panel

    # setup a GeneCollection object that contains all the genes/loci we are interested in
    wildtype_gene_collection=piezo.GeneCollection(species="M. tuberculosis",genbank_file=options.genbank,database=options.database,study=options.study,gene_panel=gene_panel)


    MUTATIONS_dict={}
    MUTATIONS_counter=0
    EFFECTS_dict={}
    EFFECTS_counter=0

    # create a copy of the wildtype_gene_collection genes which we will then alter according to the VCF file
    # need a deepcopy to ensure we take all the private variables etc with us, and point just take pointers
    sample_gene_collection=deepcopy(wildtype_gene_collection)

    (vcf_folder,vcf_filename)=os.path.split(options.vcf_file)

    # find and load the Datreant object created when the VCF was moved into place
    datreant_file=piezo.VCFMeasurement(vcf_folder)

    (n_hom,n_het,n_ref,n_null)=sample_gene_collection.apply_vcf_file(options.vcf_file)

    # (study,instance,UniqueID,site_id)=(sample_gene_collection.study,sample_gene_collection.instance,sample_gene_collection.UniqueID,sample_gene_collection.site_id)

    # store some useful features
    datreant_file.categories["GENOME_N_HOM"]=n_hom
    datreant_file.categories["GENOME_N_HET"]=n_het
    datreant_file.categories["GENOME_N_REF"]=n_ref
    datreant_file.categories["GENOME_N_NULL"]=n_null
    datreant_file.categories["SPECIES"]=sample_gene_collection.species
    datreant_file.categories["LINEAGE"]=sample_gene_collection.lineage
    datreant_file.categories["SUBLINEAGE"]=sample_gene_collection.sublineage
    datreant_file.categories["LINEAGE_PERCENTAGE"]=sample_gene_collection.lineage_percentage

    # by default assume wildtype behaviour so set all drug phenotypes to be susceptible
    phenotype={}
    for drug in walker_catalogue.drug_list:
        phenotype[drug]="S"

    # now get all the genes to calculate their own differences w.r.t the references, i.e. their mutations!
    for gene_name in tqdm(wildtype_gene_collection.gene_panel):

        # now we pass the wildtype_gene_collection gene so our VCF gene can calculate its mutations
        sample_gene_collection.gene[gene_name].identify_mutations(wildtype_gene_collection.gene[gene_name])

        mutations=sample_gene_collection.gene[gene_name].mutations

        for mutation_name in mutations:

            # store
            MUTATIONS_dict[MUTATIONS_counter]=[ vcf_filename,gene_name,\
                                                mutations[mutation_name]["VARIANT_TYPE"],\
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

                    EFFECTS_dict[EFFECTS_counter]=[vcf_filename,gene_name,mutation_name,drug_name,prediction[drug_name]]
                    EFFECTS_counter+=1
            else:
                EFFECTS_dict[EFFECTS_counter]=[vcf_filename,gene_name,mutation_name,"UNK","U"]
                EFFECTS_counter+=1

    MUTATIONS=pandas.DataFrame.from_dict(MUTATIONS_dict,orient="index",columns=["FILENAME","GENE","MUTATION_TYPE","MUTATION","ELEMENT_TYPE","POSITION","PROMOTER","CDS","SYNONYMOUS","NONSYNONYMOUS","INSERTION","DELETION","REF","ALT","NUMBER_NUCLEOTIDE_CHANGES"])
    MUTATIONS.set_index(["FILENAME","GENE","MUTATION"],inplace=True)
    MUTATIONS.to_csv(sample_gene_collection.vcf_folder+"/mutations.csv")

    EFFECTS=pandas.DataFrame.from_dict(EFFECTS_dict,orient="index",columns=["FILENAME","GENE","MUTATION","DRUG","PREDICTION"])
    EFFECTS.set_index(["FILENAME","DRUG","GENE","MUTATION"],inplace=True)
    EFFECTS.to_csv(sample_gene_collection.vcf_folder+"/effects.csv")

    wgs_prediction_string=""
    for drug in walker_catalogue.drug_list:
        datreant_file.categories["WGS_PREDICTION_"+drug]=phenotype[drug]
        wgs_prediction_string+=phenotype[drug]
    datreant_file.categories["WGS_PREDICTION_STRING"]=wgs_prediction_string

    tb_type_1="MDR"
    if datreant_file.categories["WGS_PREDICTION_INH"]=="R" and datreant_file.categories["WGS_PREDICTION_RIF"]=="R":
        if (datreant_file.categories["WGS_PREDICTION_MXF"]=="R" or datreant_file.categories["WGS_PREDICTION_LEV"]=="R") and (datreant_file.categories["WGS_PREDICTION_AMI"]=="R" or datreant_file.categories["WGS_PREDICTION_KAN"]=="R"):
            tb_type_1="XDR"
        else:
            tb_type_1="MDR"
    elif datreant_file.categories["WGS_PREDICTION_RIF"]=="R":
        tb_type_1="RIF"
    else:
        tb_type_1="SUS"

    datreant_file.categories["TB_TYPE_1"]=tb_type_1
