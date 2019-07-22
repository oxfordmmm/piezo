#! /usr/bin/python3.5

import argparse, os, pathlib
from copy import deepcopy
from datetime import datetime

import pandas
from tqdm import tqdm
from Bio import SeqIO

import piezo
import snpit
import gemucator

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_file",required=True,help="the path to a single VCF file")
    parser.add_argument("--genbank_file",default="H37Rv.gbk",help="the genbank file of the H37Rv M. tuberculosis wildtype_gene_collection genome")
    parser.add_argument("--catalogue_file",default=None,required=False,help="the path to the resistance catalogue")
    parser.add_argument("--catalogue_name",default=None,required=False,help="the name of the required catalogue, as defined in the resistance catalogue")
    parser.add_argument("--all_mutations",action='store_true',default=False,required=False,help="whether to list all mutations in a VCF file; exclusive with specifying a catalogue")
    parser.add_argument("--progress",action='store_true',default=False,help="whether to show progress using tqdm")
    parser.add_argument("--debug",action='store_true',default=False,help="print progress statements to STDOUT to help debugging")
    options = parser.parse_args()

    if options.all_mutations:
        assert (options.catalogue_file is None) and (options.catalogue_name is None), "can only specify one of --catalogue_file and --all_mutations, not both!"

    if options.catalogue_file is not None:
        assert (not options.all_mutations), "can only specify one of --catalogue_name and --all_mutations, not both!"

    tb_lineage = snpit.SnpIt(
        threshold = 10,
        ignore_filter=False,
        ignore_status=True,
    )

    results=tb_lineage.classify_vcf(options.vcf_file)

    for sample_name, (percentage, lineage) in results.items():
        output = snpit.snpit.format_output_string(sample_name, round(percentage, 2), lineage)
        print(output)

    # check the log folder exists (it probably does)
    pathlib.Path('logs/').mkdir(parents=True, exist_ok=True)

    # create a datestamp for the log files
    datestamp = datetime.strftime(datetime.now(), '%Y-%m-%d_%H%M')

    if options.catalogue_file is not None:
        if options.debug:
            print("instantiating a Resistance Catalogue..")

        # instantiate a Resistance Catalogue instance by passing a text file
        walker_catalogue=piezo.ResistanceCatalogue( input_file=options.catalogue_file,
                                                    log_file="logs/piezo-resistance-catalogue-"+datestamp+".csv",
                                                    genbank_file=options.genbank_file,
                                                    catalogue_name=options.catalogue_name )


        # retrieve the dictionary of genes from the Resistance Catalogue
        gene_panel=walker_catalogue.gene_panel

    elif options.all_mutations:

        # if no catalogue is specified, simply build up a gene_panel that is simply all gene and locus_tag
        # entries in the GenBank file
        reference=SeqIO.read(options.genbank_file,'genbank')

        gene_panel={}

        for record in reference.features:
            if record.type in ['CDS','rRNA']:
                if 'gene' in record.qualifiers.keys():
                    gene_name=record.qualifiers['gene'][0]
                    if record.type=='rRNA':
                        gene_panel[gene_name]="RNA"
                    else:
                        gene_panel[gene_name]="GENE"
                elif 'locus_tag' in record.qualifiers.keys():
                    gene_name=record.qualifiers['locus_tag'][0]
                    if record.type=='rRNA':
                        gene_panel[gene_name]="RNA"
                    else:
                        gene_panel[gene_name]="LOCUS"
                else:
                    continue

    if options.debug:
        print("creating a reference Gene Collection...")

    # setup a GeneCollection object that contains all the genes/loci we are interested in
    wildtype_gene_collection=piezo.GeneCollection(species="M. tuberculosis",genbank_file=options.genbank_file,gene_panel=gene_panel,log_file="logs/piezo-genes-"+datestamp+".csv")

    if options.debug:
        print("setting up a reference Gemucator object...")

    # setup a Gemucator object so you can check the mutations are valid
    reference_genome=gemucator.gemucator(genbank_file=options.genbank_file)

    if options.debug:
        print("creating a sample Gene Collection by copying the reference Gene Collection...")

    # create a copy of the wildtype_gene_collection genes which we will then alter according to the VCF file
    # need a deepcopy to ensure we take all the private variables etc with us, and point just take pointers
    sample_gene_collection=deepcopy(wildtype_gene_collection)

    (vcf_folder,vcf_filename)=os.path.split(options.vcf_file)

    vcf_stem=options.vcf_file.split(".vcf")[0]

    metadata={}

    if options.debug:
        print("applying the VCF file to sample Gene Collection...")

    (n_hom,n_het,n_ref,n_null)=sample_gene_collection.apply_vcf_file(options.vcf_file)

    if options.debug:
        print("working out the Lineage using snpit...")

    # call snpit to work out the specific species of Mycobacteria
    # tb=snpit.snpit(input_file=options.vcf_file)

    # store some useful features
    metadata["GENOME_N_HOM"]=n_hom
    metadata["GENOME_N_HET"]=n_het
    metadata["GENOME_N_REF"]=n_ref
    metadata["GENOME_N_NULL"]=n_null
    metadata["SNPIT_SPECIES"]=tb.species
    metadata["SNPIT_LINEAGE"]=tb.lineage
    metadata["SNPIT_SUBLINEAGE"]=tb.sublineage
    metadata["SNPIT_LINEAGE_PERCENTAGE"]="%.1f %%" % tb.percentage

    # by default assume wildtype behaviour so set all drug phenotypes to be susceptible
    phenotype={}
    for drug in walker_catalogue.drug_list:
        phenotype[drug]="S"

    MUTATIONS_dict={}
    MUTATIONS_counter=0
    EFFECTS_dict={}
    EFFECTS_counter=0

    # print(wildtype_gene_collection.gene_panel)

    # now get all the genes to calculate their own differences w.r.t the references, i.e. their mutations!
    for gene_name in wildtype_gene_collection.gene_panel:

        # now we pass the wildtype_gene_collection gene so our VCF gene can calculate its mutations
        sample_gene_collection.gene[gene_name].identify_mutations(wildtype_gene_collection.gene[gene_name])

        mutations=sample_gene_collection.gene[gene_name].mutations

        for mutation_name in mutations:

            # store
            MUTATIONS_dict[MUTATIONS_counter]=[ vcf_filename,gene_name,\
                                                mutations[mutation_name]["VARIANT_TYPE"],\
                                                mutation_name,\
                                                mutations[mutation_name]["ELEMENT_TYPE"],\
                                                int(mutations[mutation_name]["POSITION"]),\
                                                mutations[mutation_name]["PROMOTER"],\
                                                mutations[mutation_name]["CDS"],\
                                                mutations[mutation_name]["SYNONYMOUS"],\
                                                mutations[mutation_name]["NONSYNONYMOUS"],\
                                                mutations[mutation_name]["INSERTION"],\
                                                mutations[mutation_name]["DELETION"],\
                                                mutations[mutation_name]["REF"],\
                                                mutations[mutation_name]["ALT"],\
                                                mutations[mutation_name]["INDEL_1"],\
                                                mutations[mutation_name]["INDEL_2"],\
                                                mutations[mutation_name]["INDEL_3"],\
                                                int(mutations[mutation_name]["REF_COVERAGE"]),\
                                                int(mutations[mutation_name]["ALT_COVERAGE"]),\
                                                float(mutations[mutation_name]["MINOS_SCORE"]),\
                                                int(mutations[mutation_name]["NUMBER_NUCLEOTIDE_CHANGES"])]
            MUTATIONS_counter+=1

            gene_mutation=gene_name+"_"+mutation_name

            if options.catalogue_file is not None and reference_genome.valid_mutation(gene_mutation):

                prediction=walker_catalogue.predict(gene_mutation=gene_mutation,verbose=options.progress)

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


    MUTATIONS=pandas.DataFrame.from_dict(MUTATIONS_dict,orient="index",columns=["FILENAME","GENE","MUTATION_TYPE","MUTATION","ELEMENT_TYPE","POSITION","PROMOTER","CDS","SYNONYMOUS","NONSYNONYMOUS","INSERTION","DELETION","REF","ALT","INDEL_1","INDEL_2","INDEL_3","REF_COVERAGE","ALT_COVERAGE","MINOS_SCORE","NUMBER_NUCLEOTIDE_CHANGES"])
    MUTATIONS.set_index(["FILENAME","GENE","MUTATION"],inplace=True)
    MUTATIONS.to_csv(vcf_stem+"-mutations.csv",header=True)

    if options.catalogue_file is not None:
        EFFECTS=pandas.DataFrame.from_dict(EFFECTS_dict,orient="index",columns=["FILENAME","GENE","MUTATION","DRUG","PREDICTION"])
        EFFECTS.set_index(["FILENAME","DRUG","GENE","MUTATION"],inplace=True)
        EFFECTS.to_csv(vcf_stem+"-effects.csv",header=True)

        wgs_prediction_string=""
        for drug in walker_catalogue.drug_list:
            metadata["WGS_PREDICTION_"+drug]=phenotype[drug]

    print("%40s %s" % ("VCF file", options.vcf_file))
    if options.catalogue_file is not None:
        print("%40s %s" % ("Catalogue", options.catalogue_file))
        print("%40s %s" % ("Catalogue_name", options.catalogue_name))
    print("%40s %s" % ("Genbank file", options.genbank_file))
    for i in sorted(metadata):
        print("%40s %s" % (i, metadata[i]))
