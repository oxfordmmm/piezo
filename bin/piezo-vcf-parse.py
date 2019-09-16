#! /usr/bin/python3.5

import argparse, os, pathlib, pickle, gzip

from copy import deepcopy
from datetime import datetime
from collections import defaultdict

import pandas, numpy
from tqdm import tqdm

import gumpy, piezo

def variant_infer_columns(row):
    if row["IS_INDEL"]:
        ref=sample_genome.indel_ref[sample_genome.genome_index==row["INDEX"]][0]
        alt=sample_genome.indel_alt[sample_genome.genome_index==row["INDEX"]][0]
        indel_length=sample_genome.indel_length[sample_genome.genome_index==row["INDEX"]][0]
        variant=str(row['INDEX'])+"_indel"
    else:
        ref=reference_genome.genome_sequence[sample_genome.genome_index==row["INDEX"]][0]
        alt=sample_genome.genome_sequence[sample_genome.genome_index==row["INDEX"]][0]
        indel_length=0
        variant=str(row['INDEX'])+ref+">"+alt
    return pandas.Series([ref,alt,indel_length,variant])

def variant_assign_booleans(row):
    in_promoter=False
    in_cds=False
    associated_with_gene=False
    is_insertion=False
    is_deletion=False
    if row["POSITION"]<0:
        in_promoter=True
    if row["POSITION"]>0:
        in_cds=True
    if row["POSITION"]!=0:
        associated_with_gene=True
    return(pandas.Series([in_promoter,in_cds,associated_with_gene]))

def mutations_assign_booleans(row):
    if row['POSITION']<0:
        is_cds=False
        is_promoter=True
    else:
        is_cds=True
        is_promoter=False
    if row['MUTATION'][0].isupper():
        if row['MUTATION'][0]==row['MUTATION'][-1]:
            is_synonymous=True
            is_nonsynonymous=False
        else:
            is_synonymous=False
            is_nonsynonymous=True
    else:
        is_synonymous=False
        is_nonsynonymous=False

    is_het=False
    is_null=False
    is_snp=False

    if row["MUTATION"][-1] in ['z','Z']:
        is_het=True
    elif row["MUTATION"][-1] in ['x','X']:
        is_null=True
    elif not row["IS_INDEL"]:
        is_snp=True

    return(pandas.Series([is_promoter,is_cds,is_synonymous,is_nonsynonymous,is_het,is_null,is_snp]))

def mutations_count_number_nucleotide_changes(row):
    if row['REF'] is not None and len(row['REF'])==3:
        return(numpy.sum([i!=j for (i,j) in zip(row['REF'],row['ALT'] ) ]))
    else:
        return(0)

def mutations_split_mutation(row):
    if "indel" in row["MUTATION"]:
        cols=row["MUTATION"].split("_")
        return(None,int(cols[0]),None,False,True,int(cols[2]),"INDEL")
    else:
        mut=row["MUTATION"]

        if mut[0].isupper():
            gene=row["GENE"]
            ref=reference_genome.genes[gene].codons[sample_genome.genes[gene].amino_acid_numbering==int(mut[1:-1])][0]
            alt=sample_genome.genes[gene].codons[sample_genome.genes[gene].amino_acid_numbering==int(mut[1:-1])][0]
            return pandas.Series([ref,int(mut[1:-1]),alt,True,False,None,"SNP"])
        else:
            return pandas.Series([mut[0],int(mut[1:-1]),mut[-1],True,False,None,"SNP"])

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--vcf_file",required=True,help="the path to a single VCF file")
    parser.add_argument("--genome_object",default="H37Rv_3.pkl.gz",help="the path to a compressed pickled gumpy Genome object")
    parser.add_argument("--catalogue_file",default=None,required=False,help="the path to the resistance catalogue")
    parser.add_argument("--catalogue_name",default=None,required=False,help="the name of the required catalogue, as defined in the resistance catalogue")
    parser.add_argument("--ignore_vcf_status",action='store_true',default=False,help="whether to ignore the STATUS field in the vcf (e.g. necessary for some versions of Clockwork VCFs)")
    parser.add_argument("--ignore_vcf_filter",action='store_true',default=False,help="whether to ignore the FILTER field in the vcf (e.g. necessary for some versions of Clockwork VCFs)")
    parser.add_argument("--progress",action='store_true',default=False,help="whether to show progress using tqdm")
    parser.add_argument("--debug",action='store_true',default=False,help="print progress statements to STDOUT to help debugging")
    options = parser.parse_args()

    if options.catalogue_file:
        assert options.catalogue_name is not None, "if you specify a catalogue_file you must also specify a catalogue name!"

    if options.debug:
        print("Loading the reference object..")

    # load the pickled reference gumpy Genome object to save time
    INPUT=gzip.open(options.genome_object,"rb")
    reference_genome=pickle.load(INPUT)
    INPUT.close()

    # check the log folder exists (it probably does)
    pathlib.Path('logs/').mkdir(parents=True, exist_ok=True)

    # create a datestamp for the log files
    datestamp = datetime.strftime(datetime.now(), '%Y-%m-%d_%H%M')

    if options.catalogue_file is not None:
        if options.debug:
            print("Instantiating a Resistance Catalogue..")

        # instantiate a Resistance Catalogue instance by passing a text file
        resistance_catalogue=piezo.ResistanceCatalogue( input_file=options.catalogue_file,
                                                    log_file="logs/piezo-resistance-catalogue-"+datestamp+".csv",
                                                    gumpy_genome=reference_genome,
                                                    catalogue_name=options.catalogue_name )

    if options.debug:
        print("Creating a sample Genome object by copying the reference Genome object...")

    # create a copy of the reference genome which we will then alter according to the VCF file
    # need a deepcopy to ensure we take all the private variables etc with us
    sample_genome=deepcopy(reference_genome)

    (vcf_folder,vcf_filename)=os.path.split(options.vcf_file)

    vcf_stem=options.vcf_file.split(".vcf")[0]

    metadata={}

    if options.debug:
        print("applying the VCF file to sample Genome...")

    # FIXME:(n_hom,n_het,n_ref,n_null)=
    sample_genome.apply_vcf_file( show_progress_bar=options.progress,\
                                  vcf_file=options.vcf_file,\
                                  ignore_status=True,\
                                  ignore_filter=False,\
                                  metadata_fields=['GT_CONF','GT_CONF_PERCENTILE'])

    # #FIXME store some useful features
    # metadata["GENOME_N_HOM"]=n_hom
    # metadata["GENOME_N_HET"]=n_het
    # metadata["GENOME_N_REF"]=n_ref
    # metadata["GENOME_N_NULL"]=n_null


    if options.debug:
        print("Creating the VARIANTS table..")

    # create a dictionary
    VARIANT_dict=defaultdict(str)

    # return the list of genome indices where there is variation (SNPs and INDELs)
    index=reference_genome-sample_genome

    # create a Boolean mask for fancy indexing
    mask=numpy.isin(sample_genome.genome_index,index)

    # create the columns that are common to both SNPs and INDELs for all variants
    VARIANT_dict['INDEX']=index
    VARIANT_dict['GENE']=sample_genome.genome_feature_name[mask]
    VARIANT_dict['POSITION']=sample_genome.genome_positions[mask]
    VARIANT_dict['NUMBERING']=sample_genome.genome_numbering[mask]
    VARIANT_dict['IS_SNP']=sample_genome.is_snp[mask]
    VARIANT_dict['IS_HET']=sample_genome.is_het[mask]
    VARIANT_dict['IS_INDEL']=sample_genome.is_indel[mask]
    VARIANT_dict['IS_NULL']=sample_genome.is_null[mask]
    VARIANT_dict['COVERAGE']=sample_genome.coverage[mask]
    VARIANT_dict['GT_CONF']=sample_genome.genome_sequence_metadata["GT_CONF"][mask]
    VARIANT_dict['GT_CONF_PERCENTILE']=sample_genome.genome_sequence_metadata["GT_CONF_PERCENTILE"][mask]
    VARIANT_dict["HET_VARIANT_0"]=sample_genome.het_variations[mask][:,0]
    VARIANT_dict["HET_VARIANT_1"]=sample_genome.het_variations[mask][:,1]
    VARIANT_dict["HET_COVERAGE_0"]=sample_genome.het_coverage[mask][:,0]
    VARIANT_dict["HET_COVERAGE_1"]=sample_genome.het_coverage[mask][:,1]
    VARIANT_dict["HET_INDEL_LENGTH_0"]=sample_genome.het_indel_length[mask][:,0]
    VARIANT_dict["HET_INDEL_LENGTH_1"]=sample_genome.het_indel_length[mask][:,1]
    VARIANT_dict["HET_REF"]=sample_genome.het_ref[mask]
    VARIANT_dict["HET_ALT_0"]=sample_genome.het_alt[mask][:,0]
    VARIANT_dict["HET_ALT_1"]=sample_genome.het_alt[mask][:,1]

    # create a preliminary dataframe
    VARIANT=pandas.DataFrame(data=VARIANT_dict)

    # populate the key columns
    VARIANT[['REF','ALT','INDEL_LENGTH','VARIANT']]=VARIANT.apply(variant_infer_columns,axis=1)

    VARIANT['UNIQUEID']=sample_genome.name

    # define some Boolean columns for ease of analysis
    VARIANT[["IN_PROMOTER","IN_CDS","ASSOCIATED_WITH_GENE"]]=VARIANT.apply(variant_assign_booleans,axis=1)

    # reorder the columns
    VARIANT=VARIANT[['UNIQUEID','INDEX','VARIANT','ASSOCIATED_WITH_GENE','GENE','NUMBERING','POSITION','IS_HET','IS_INDEL','IS_NULL','IS_SNP','INDEL_LENGTH','IN_PROMOTER','IN_CDS','REF','ALT','COVERAGE','GT_CONF','GT_CONF_PERCENTILE','HET_ALT_0','HET_ALT_1','HET_COVERAGE_0','HET_COVERAGE_1','HET_INDEL_LENGTH_0','HET_INDEL_LENGTH_1','HET_REF','HET_VARIANT_0','HET_VARIANT_1']]

    # set the index
    VARIANT.set_index(['UNIQUEID','VARIANT'],inplace=True,verify_integrity=True)

    # save to a CSV file
    VARIANT.to_csv(vcf_stem+"-VARIANTS.csv",header=True)

    if options.debug:
        print("Creating the MUTATIONS table..")

    # create a dictionary
    MUTATIONS_dict=defaultdict(str)

    # these are the only three columns we can populate quickly
    for i in ['GENE','MUTATION','ELEMENT_TYPE']:
        MUTATIONS_dict[i]=[]

    for gene in reference_genome.gene_names:

        # find out the mutations
        mutations=sample_genome.genes[gene].list_mutations_wrt(reference_genome.genes[gene])

        if mutations is not None:

            number=len(mutations)

            MUTATIONS_dict["GENE"]=numpy.append(MUTATIONS_dict["GENE"],number*[gene])
            MUTATIONS_dict["MUTATION"]=numpy.append(MUTATIONS_dict["MUTATION"],mutations)
            gene_type=reference_genome.genes[gene].gene_type
            MUTATIONS_dict['ELEMENT_TYPE']=numpy.append(MUTATIONS_dict['ELEMENT_TYPE'],number*[gene_type])

    # create a preliminary dataframe
    MUTATIONS=pandas.DataFrame(data=MUTATIONS_dict)

    # infer the remaining key columns
    MUTATIONS[['REF','POSITION','ALT','IS_SNP','IS_INDEL','INDEL_LENGTH',"MUTATION_TYPE"]] = MUTATIONS.apply(mutations_split_mutation,axis=1)

    # define some Boolean columns for ease of analysis
    MUTATIONS[['IN_PROMOTER','IN_CDS','IS_SYNONYMOUS','IS_NONSYNONYMOUS','IS_HET','IS_NULL','IS_SNP']]=MUTATIONS.apply(mutations_assign_booleans,axis=1)

    # calculate the number of nucleotide changes required for this mutation
    MUTATIONS["NUMBER_NUCLEOTIDE_CHANGES"]=MUTATIONS.apply(mutations_count_number_nucleotide_changes,axis=1)

    MUTATIONS['UNIQUEID']=sample_genome.name

    # reorder the columns
    MUTATIONS=MUTATIONS[["UNIQUEID","GENE","MUTATION","MUTATION_TYPE","ELEMENT_TYPE","POSITION","IN_PROMOTER","IN_CDS","IS_SYNONYMOUS","IS_NONSYNONYMOUS","IS_HET","IS_INDEL","IS_NULL","IS_SNP","REF","ALT","NUMBER_NUCLEOTIDE_CHANGES"]]

    # set the index
    MUTATIONS.set_index(["UNIQUEID","GENE",'MUTATION'],inplace=True,verify_integrity=True)

    # save to a CSV file
    MUTATIONS.to_csv(vcf_stem+"-MUTATIONS.csv")

    MUTATIONS.reset_index(inplace=True)

    # can only infer predicted effects and ultimate phenotypes if a resistance catalogue has been supplied!
    if options.catalogue_file is not None:

        # subset down to only those mutations in the catalogue for making predictions
        MUTATIONS_IN_CATALOGUE=MUTATIONS.loc[MUTATIONS.GENE.isin(resistance_catalogue.gene_list)]

        # by default assume wildtype behaviour so set all drug phenotypes to be susceptible
        phenotype={}
        for drug in resistance_catalogue.drug_list:
            phenotype[drug]="S"

        EFFECTS_dict={}
        EFFECTS_counter=0

        for gene_name,mutation_name in zip(MUTATIONS_IN_CATALOGUE.GENE,MUTATIONS_IN_CATALOGUE.MUTATION):

            try:
                valid=reference_genome.valid_gene_mutation(gene_name+"_"+mutation_name)
            except:
                valid=False

            assert valid, gene_name+"_"+mutation_name+" is not a valid mutation!"

            prediction=resistance_catalogue.predict(gene_mutation=gene_name+"_"+mutation_name,verbose=False)

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

                    EFFECTS_dict[EFFECTS_counter]=[sample_genome.name,gene_name,mutation_name,options.catalogue_name,drug_name,prediction[drug_name]]
                    EFFECTS_counter+=1
            else:
                EFFECTS_dict[EFFECTS_counter]=[sample_genome.name,gene_name,mutation_name,options.catalogue_name,"UNK","S"]
                EFFECTS_counter+=1

        EFFECTS=pandas.DataFrame.from_dict(EFFECTS_dict,orient="index",columns=["UNIQUEID","GENE","MUTATION","CATALOGUE_NAME","DRUG","PREDICTION"])
        EFFECTS.set_index(["UNIQUEID","DRUG","GENE","MUTATION","CATALOGUE_NAME"],inplace=True)
        EFFECTS.to_csv(vcf_stem+"-EFFECTS.csv")

        wgs_prediction_string=""
        for drug in resistance_catalogue.drug_list:
            metadata["WGS_PREDICTION_"+drug]=phenotype[drug]

    print("%40s %s" % ("VCF file", options.vcf_file))
    if options.catalogue_file is not None:
        print("%40s %s" % ("Catalogue", options.catalogue_file))
        print("%40s %s" % ("Catalogue_name", options.catalogue_name))
    print("%40s %s" % ("Genome file", options.genome_object))
    for i in sorted(metadata):
        print("%40s %s" % (i, metadata[i]))
