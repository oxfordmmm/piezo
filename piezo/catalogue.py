#! /usr/bin/env python

import logging, pkg_resources, pathlib
from datetime import datetime

import pandas, numpy

from piezo import Gene

class ResistanceCatalogue(object):

    def __init__(self,input_file=None):

        # read in the Walker Resistance Catalogue and make a dictionary
        self.entry={}

        # define the list of drugs
        # ordered as
        # 1st line: R,I,P,E
        # 2nd line: aminoglycosides (AMI, KAN), fluroquinolones (LEV, MXF), ETH, PAS
        # 3rd line: RFB, LZD, BDQ, DLM
        # Repurposed: CFZ
        self.drug_list=["RIF","INH","PZA","EMB","AMI","KAN","LEV","MXF","ETH","PAS","RFB","LZD","BDQ","DLM","CFZ"]

        # remember the config path
        self.config_path = '/'.join(('..','config'))

        self._parse_catalogue_file(pkg_resources.resource_filename("cryptic", self.config_path+"/"+input_file))

        # check the log folder exists (it probably does)
        pathlib.Path('logs/').mkdir(parents=True, exist_ok=True)


        datestamp = datetime.strftime(datetime.now(), '%Y-%m-%d_%H%M')
        logging.basicConfig(filename="logs/cryptic-genetics-resistancecatalogue-"+datestamp+".csv",level=logging.INFO,format='%(levelname)s, %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

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
        # print((self.gene_panel))

        # print(self.gene_list)

        # self.resistance_catalogue.set_index(["DRUG","MUTATION"],inplace=True)


    # def _parse_catalogue_file(self,input_file):
    #
    #     self.gene_lookup={}
    #     self.drug_lookup={}
    #
    #     # open the catalogue file
    #     with open(input_file) as file:
    #
    #         # iterate through it line-by-line
    #         for line in file:
    #
    #             # split on commas
    #             cols=line.split(',')
    #
    #             # skip the header row
    #             if cols[0]!="DRUG":
    #
    #                 # if marked as a PROTein mutation
    #                 # FIXME cannot handle rrs RNA mutations
    #                 assert cols[3] in ['CDS','PROM','RNA'], cols[3]+' promoter location other than CDS/PROM/RNA specified in the catalogue file!'
    #
    #                 # FIXME cannot handle indels
    #                 if cols[2]=="SNP" and "Stop" not in cols[1]:
    #
    #                     assert cols[4] in ['R','S','U'], 'phenotype other than R/S/U specified in the catalogue file!'
    #
    #                     drug_name=cols[0]
    #                     gene_name=cols[1].split("_")[0]
    #
    #                     mutation=cols[1].split("_")[1]
    #
    #                     position=int(mutation[1:-1])
    #
    #                     if position<0 or cols[3]=='RNA':
    #                         before=mutation[:1].lower()
    #                         after=mutation[-1:].lower()
    #                     else:
    #                         before=mutation[:1].upper()
    #                         after=mutation[-1:].upper()
    #
    #                     phenotype=cols[4]
    #
    #                     if drug_name not in self.drug_lookup.keys():
    #                         self.drug_lookup[drug_name]=[]
    #                     else:
    #                         if gene_name not in self.drug_lookup[drug_name]:
    #                             self.drug_lookup[drug_name].append(gene_name)
    #
    #                     if gene_name not in self.gene_lookup.keys():
    #                         self.gene_lookup[gene_name]=[]
    #                     else:
    #                         if drug_name not in self.gene_lookup[gene_name]:
    #                             self.gene_lookup[gene_name].append(drug_name)
    #
    #                     if gene_name not in self.entry.keys():
    #                         self.entry[gene_name]={}
    #                     if position not in self.entry[gene_name].keys():
    #                         self.entry[gene_name][position]={}
    #                     if drug_name not in self.entry[gene_name][position].keys():
    #                         self.entry[gene_name][position][drug_name]={}
    #
    #                     if 'after' in self.entry[gene_name][position][drug_name].keys():
    #
    #                         # check that the amino acid at this position according to the genbank file matches what is in the catalogue
    #                         if before!=self.entry[gene_name][position][drug_name]['before']:
    #                             logging.warning("before is "+before+" at position "+str(position)+" in gene "+gene_name+", not "+self.entry[gene_name][position][drug_name]['before'])
    #
    #                         # check this amino acid doesn't already exist!
    #                         if after in self.entry[gene_name][position][drug_name]['after']:
    #                             logging.warning("duplicate row in catalogue: "+line.rstrip()) #,gene_name,position,after,self.entry[gene_name][position][after],self.entry[gene_name][position][after].keys())
    #
    #                         # if it doesn't, append it to the lists
    #                         else:
    #                             self.entry[gene_name][position][drug_name]['after'].append(after)
    #                             self.entry[gene_name][position][drug_name]['phenotype'].append(phenotype)
    #
    #                     else:
    #                         self.entry[gene_name][position][drug_name]['before']=before
    #                         self.entry[gene_name][position][drug_name]['after']=[after]
    #                         self.entry[gene_name][position][drug_name]['phenotype']=[phenotype]
    #                         # self.entry[gene_name][position][drug_name]={'before':before,'phenotype':phenotype,'after':after}
    #                     # else:
    #                     #     print("duplicate row in catalogue: "+line.rstrip()) #,gene_name,position,after,self.entry[gene_name][position][after],self.entry[gene_name][position][after].keys())
    #
    def predict(self,mutation=None,catalogue_name="CRYPTICv0.9"):

        assert catalogue_name+"_PREDICTION" in self.resistance_catalogue.columns, "specified catalogue does not have a column in the loaded Resistance Catalogue!"

        components=mutation.split("_")

        gene_name=components[0]

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
            assert components[2] in ['ins','del'], "INDEL can only be ins or del, given "+mutation
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
                        matching_rows = self.resistance_catalogue.loc[(self.resistance_catalogue.GENE_MUTATION==gene_name+"_"+before+str(position)+"*") & (self.resistance_catalogue.DRUG==compound)]
                        if len(matching_rows)==1:
                            result[compound]=matching_rows[catalogue_name+"_PREDICTION"].values[0]
                            found_result=True
                    elif mutation_type=="INDEL":
                        matching_rows = self.resistance_catalogue.loc[(self.resistance_catalogue.GENE_MUTATION==gene_name+"_"+str(position)+"_"+indel_type+"_*") & (self.resistance_catalogue.DRUG==compound)]
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

# def predict(self,gene_name=None,before=None,position=None,after=None):
#
#     result={}
#
#     # if a promoter mutation
#     if position<0:
#
#         # check isn't lower than -100
#         assert position>=-100, "promoter mutation too far from CDS; position is "+str(position)
#
#         # and insist it is a nucleotide
#         assert after in ['a','c','t','g','x','z'], after+" is not a base! "+gene_name+" "+str(position)
#
#     # otherwise it is a mutation in the coding region
#     else:
#
#         # insist we've been given an amino acid or a wildcard only
#         assert after in ['a','c','t','g','x','z','*','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid!"
#
#     # if gene is not in our panel, then return an S as we just don't know
#     if gene_name not in self.gene_list:
#
#         # if it isn't, then we don't know what drugs it might affect, so just return susceptible
#         return("S")
#
#     # otherwise it _must_ be in our catalogue
#     else:
#
#         # find out the drugs to which changes in this gene can confer resistance
#         drugs=self.gene_lookup[gene_name]
#
#         # if there aren't any just return an S again
#         if drugs==[]:
#             return("S")
#
#         # create the result dictionary e.g. {"RIF":'R', "RFB":'R'}
#         # most of the time the predictions will be identical, but this allows them to diverge
#         result={}
#
#         # deal with each compound, one at a time
#         for compound in drugs:
#
#             # if the mutation is synonmous, just return SUSCEPTIBLE
#             if before==after:
#                 result[compound]="S"
#
#             # if there are no recorded mutations at that position for any drugs, return UNKNOWN
#             elif position not in self.entry[gene_name].keys():
#                 result[compound]="U"
#
#             # otherwise there must be at least one mutation in that gene at this position
#             elif position in self.entry[gene_name].keys():
#
#                 # first check that there is an entry for this drug
#                 if compound in self.entry[gene_name][position].keys():
#
#                     # if there is, first check the stored before and passed before match
#                     if self.entry[gene_name][position][compound]['before']!=before:
#                         logging.warning("Mismatch between catalogue and genbank file for wildtype amino acid in "+gene_name+" at position "+str(position)+" genbank:"+before+" catalog:"+self.entry[gene_name][position][compound]['before'])
#
#                     found_result=False
#
#                     # if there is, check if either a wildcard or exact match is specified
#                     if "*" in self.entry[gene_name][position][compound]['after']:
#
#                         # then find out the position in the list
#                         index=self.entry[gene_name][position][compound]['after'].index("*")
#
#                         # then assign the stored phenotype from the correct position in the list
#                         result[compound]=self.entry[gene_name][position][compound]['phenotype'][index]
#
#                         found_result=True
#
#                     # otherwise look for an exact match
#                     if after in self.entry[gene_name][position][compound]['after']:
#
#                         # then find out the position in the list
#                         index=self.entry[gene_name][position][compound]['after'].index(after)
#
#                         # then assign the stored phenotype from the correct position in the list
#                         result[compound]=self.entry[gene_name][position][compound]['phenotype'][index]
#
#                         found_result=True
#
#                     # otherwise there is no match in the catalogue at this position, return UNKNOWN
#                     if not found_result:
#                         result[compound]="U"
#
#                 # otherwise there is no mutation for this drug at this position, return UNKNOWN
#                 else:
#                     result[compound]="U"
#
#         return(result)
