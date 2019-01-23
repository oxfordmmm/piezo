#! /usr/bin/env python

import logging, pkg_resources

import pandas, numpy

from piezo import Gene
from gemucator import gemucator

class ResistanceCatalogue(object):

    def __init__(self,input_file=None,log_file=None,genbank_file=None):

        # read in the Walker Resistance Catalogue and make a dictionary
        self.entry={}

        # remember the config path
        self.config_path = '/'.join(('..','config'))

        self.reference_genome=gemucator(genbank_file=genbank_file)

        # self._parse_catalogue_file(pkg_resources.resource_filename("piezo", self.config_path+"/"+input_file))
        self._parse_catalogue_file(input_file)

        logging.basicConfig(filename=log_file,level=logging.INFO,format='%(levelname)s, %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

    def _parse_catalogue_file(self,input_file):

        self.gene_lookup={}
        self.drug_lookup={}

        self.resistance_catalogue = pandas.read_csv(input_file)

        self.number_rows = len(self.resistance_catalogue)

        # be defensive
        assert numpy.sum(self.resistance_catalogue["VARIANT_TYPE"].isin([numpy.nan,"SNP","INDEL"]))==self.number_rows, "TYPE column contains entries other than SNP or INDEL"
        assert numpy.sum(self.resistance_catalogue["VARIANT_AFFECTS"].isin([numpy.nan,'CDS','PROM','RNA']))==self.number_rows, "AFFECTS column contains entries other than CDS, PROM or RNA"
        # assert numpy.sum(self.resistance_catalogue["PRED"].isin(['R','S','U']))==self.number_rows, "PRED column contains entries other than R, S or U"

        n_duplicated_rows=len(self.resistance_catalogue.loc[self.resistance_catalogue.duplicated(subset=["DRUG","GENE","MUTATION"],keep='first')])
        assert n_duplicated_rows==0, "There are duplicated rows in the catalogue!"

        # self.resistance_catalogue=self.resistance_catalogue.loc[~self.resistance_catalogue.duplicated(subset=["DRUG","GENE","MUTATION"],keep='first')]

        self.resistance_catalogue["GENE_MUTATION"]=self.resistance_catalogue["GENE"]+"_"+self.resistance_catalogue["MUTATION"]

        for gene_name in self.resistance_catalogue.GENE.unique():
            tmp=self.resistance_catalogue.loc[(self.resistance_catalogue.GENE==gene_name) & (~self.resistance_catalogue.DRUG.isna())]
            if tmp.empty:
                self.gene_lookup[gene_name]=None
            else:
                self.gene_lookup[gene_name]=list(tmp.DRUG.unique())

        for drug_name in self.resistance_catalogue.DRUG.unique():
            tmp=self.resistance_catalogue.loc[self.resistance_catalogue.DRUG==drug_name]
            self.drug_lookup[drug_name]=list(tmp.GENE.unique())

        self.gene_list=list(self.resistance_catalogue.GENE.unique())

        self.gene_panel={}
        foo=self.resistance_catalogue[["GENE","GENE_TYPE","MUTATION"]].groupby(["GENE","GENE_TYPE"]).count().to_dict()
        for i in foo['MUTATION']:
            self.gene_panel[i[0]]=i[1]


    def predict(self,mutation=None,catalogue_name="CRYPTICv0.9"):

        assert catalogue_name+"_PREDICTION" in self.resistance_catalogue.columns, "specified catalogue does not have a column in the loaded Resistance Catalogue!"

        components=mutation.split("_")

        gene_name=components[0]

        # check that the gene exists!
        assert reference_genome.valid_gene(gene_name), gene_name+" does not exist in the specified GENBANK file!"

        # check that the reference is valid
        assert reference_genome.valid_mutation(mutation), "gene exists but "+mutation+" is badly formed; check the reference amino acid or nucleotide!"

        if len(components)==2:
            mutation_type="SNP"
            before=components[1][:1]
            after=components[1][-1:]
            position=int(components[1][1:-1])

            # insist we've been given an amino acid or a wildcard only
            assert after in ['a','c','t','g','x','z','*',"!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid or base!"
            assert before in ['a','c','t','g','x','z','*','!','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], before+" is not an amino acid or base!"

            # check the position isn't lower than -100
            assert position>=-100, "promoter mutation too far from CDS; max allowed is -100, position is "+str(position)


        elif len(components)==4:
            mutation_type="INDEL"
            position=int(components[1])
            assert components[2] in ['ins','del','indel','fs'], "INDEL can only be indel,ins,del or fs given "+mutation
            indel_type=components[2]
        else:
            raise ValueError("malformed mutation passed to predict: "+mutation)

        if gene_name not in self.gene_list:

            # if it isn't, then we don't know what drugs it might affect, so just return susceptible
            return("S")

        # otherwise it _must_ be in our catalogue
        else:

            # find out the drugs to which changes in this gene can confer resistance
            drugs=self.gene_lookup[gene_name]

            # if there aren't any just return an S again
            if drugs is None:
                return("S")

            # create the result dictionary e.g. {"RIF":'R', "RFB":'R'}
            # most of the time the predictions will be identical, but this allows them to diverge
            result={}

            # deal with each compound, one at a time
            for compound in drugs:

                same_position_rows = self.resistance_catalogue.loc[(self.resistance_catalogue.POSITION==position) & (self.resistance_catalogue.DRUG==compound)]

                # if the mutation is synonmous, just return SUSCEPTIBLE
                if mutation_type=="SNP" and before==after:
                    result[compound]="S"

                elif len(same_position_rows)==0:
                    result[compound]="U"

                else:
                    found_result=False

                    # check for wildtype matching rows
                    if mutation_type=="SNP":
                        matching_rows = self.resistance_catalogue.loc[(self.resistance_catalogue.GENE_MUTATION==gene_name+"_"+before+str(position)+"?") & (self.resistance_catalogue.DRUG==compound)]
                        if len(matching_rows)==1:
                            result[compound]=matching_rows[catalogue_name+"_PREDICTION"].values[0]
                            found_result=True
                    elif mutation_type=="INDEL":
                        matching_rows = self.resistance_catalogue.loc[(self.resistance_catalogue.GENE_MUTATION==gene_name+"_"+str(position)+"_"+indel_type+"_?") & (self.resistance_catalogue.DRUG==compound)]
                        if len(matching_rows)==1:
                            result[compound]=matching_rows[catalogue_name+"_PREDICTION"].values[0]
                            found_result=True

                    # otherwise look for an exact match
                    matching_rows = self.resistance_catalogue.loc[(self.resistance_catalogue.GENE_MUTATION==mutation) & (self.resistance_catalogue.DRUG==compound)]
                    if len(matching_rows)==1:
                        result[compound]=matching_rows[catalogue_name+"_PREDICTION"].values[0]
                        found_result=True

                    # otherwise there is no match in the catalogue at this position, return UNKNOWN
                    if not found_result:
                        result[compound]="U"

        return(result)
