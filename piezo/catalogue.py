#! /usr/bin/env python

import logging, os

import pandas, numpy

from piezo import Gene
from gemucator import gemucator

class ResistanceCatalogue(object):

    def __init__(self,input_file=None,log_file=None,genbank_file=None,catalogue_name="LID2015B"):

        '''
        Instantiate a ResistanceCatalogue

        Args:
            input_file (str): path to a resistance catalogue as a CSV file
            genbank_file (str): path to the matching GenBank file describing the reference genome
            log_file (str): path to a logfile
        '''

        self.catalogue_name=catalogue_name

        # read in the Walker Resistance Catalogue and make a dictionary
        self.entry={}

        # instantiate a gemucator instance using the same GenBank file so we can validate the mutations later on
        self.reference_genome=gemucator(genbank_file=os.path.abspath(genbank_file))

        # read in the catalogue file
        self._parse_catalogue_file(input_file)

        # start logging using the provided log_file path
        logging.basicConfig(filename=log_file,level=logging.INFO,format='%(levelname)s, %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

    def _parse_catalogue_file(self,input_file):

        '''
        Read in the Antimicrobial Resistance Catalogue.

        Args:
            input_file (str): path to a resistance catalogue as a CSV file

        Notes:
            * Applies some checks to ensure the catalogue is in the right format.
        '''

        # create some empty dictionaries
        self.gene_lookup={}
        self.drug_lookup={}

        # use Pandas to read the catalogue into a DataFrame
        self.resistance_catalogue = pandas.read_csv(input_file)

        # ensure that the specified catalogue has a column
        assert self.catalogue_name+"_PREDICTION" in self.resistance_catalogue.columns, "specified catalogue "+self.catalogue_name+" does not have a column in the loaded Resistance Catalogue!"

        assert len(self.resistance_catalogue.GENBANK_REFERENCE.unique())==1, "multiple genbank references specified in catalogue!"

        catalogue_version=self.resistance_catalogue.GENBANK_REFERENCE.unique()[0]

        assert catalogue_version==self.reference_genome.version, "Catalogue uses version "+catalogue_version+" whilst genbank file is version "+self.reference_genome.version+" !!"

        # only pull out genes which have at least one R row
        relevant_genes=self.resistance_catalogue.loc[self.resistance_catalogue[self.catalogue_name+"_PREDICTION"].isin(['S','U','R'])].GENE.unique()

        # downsample the catalogue to include include rows which have definite R/S/U and genes which have at least 1 R row
        self.resistance_catalogue=self.resistance_catalogue.loc[(self.resistance_catalogue[self.catalogue_name+"_PREDICTION"].notna()) & (self.resistance_catalogue.GENE.isin(relevant_genes))]

        # find out the number of rows
        self.number_rows = len(self.resistance_catalogue)

        # be defensive and check these columns only contain what we think they should contain
        assert numpy.sum(self.resistance_catalogue["VARIANT_TYPE"].isin([numpy.nan,"SNP","INDEL"]))==self.number_rows, "TYPE column contains entries other than SNP or INDEL"
        assert numpy.sum(self.resistance_catalogue["VARIANT_AFFECTS"].isin([numpy.nan,'CDS','PROM','RNA']))==self.number_rows, "AFFECTS column contains entries other than CDS, PROM or RNA"

        # insist that there are no duplicated rows
        n_duplicated_rows=len(self.resistance_catalogue.loc[self.resistance_catalogue.duplicated(subset=["DRUG","GENE","MUTATION"],keep='first')])
        assert n_duplicated_rows==0, "There are duplicated rows in the catalogue!"

        # iterate through the gene names
        for gene_name in self.resistance_catalogue.GENE.unique():

            # be defensive and check this gene exists in our reference
            assert self.reference_genome.valid_gene(gene_name), gene_name+" does not exist in the supplied GenBank file!"

            # select the rows matching this gene where a drug is specified
            tmp=self.resistance_catalogue.loc[(self.resistance_catalogue.GENE==gene_name) & (~self.resistance_catalogue.DRUG.isna())]

            # if there are none, then this gene has been added to the catalogue without any associated drugs as we merely want to track its mutations
            if tmp.empty:

                # in which case, record that no drugs are associated
                self.gene_lookup[gene_name]=None

            else:

                # otherwise record the drugs whose action can be affected by variation in this gene
                self.gene_lookup[gene_name]=list(tmp.DRUG.unique())

        # now we do it the other way round!
        for drug_name in self.resistance_catalogue.DRUG.unique():

            # select the rows which match this drug
            tmp=self.resistance_catalogue.loc[self.resistance_catalogue.DRUG==drug_name]

            # find out and record the genes
            self.drug_lookup[drug_name]=list(tmp.GENE.unique())

        # record all the genes present in the catalogue as a simple list
        # (this can include genes which have 'blank' rows for which no drug is yet associated)
        self.gene_list=list(self.resistance_catalogue.GENE.unique())

        self.drug_list=list(self.resistance_catalogue.DRUG.unique())

        # finally, make the gene_panel dictionary which tells you whether the gene is a LOCUS or a GENE in the reference
        self.gene_panel={}
        foo=self.resistance_catalogue[["GENE","GENE_TYPE","MUTATION"]].groupby(["GENE","GENE_TYPE"]).count().to_dict()
        for i in foo['MUTATION']:
            self.gene_panel[i[0]]=i[1]

    def predict(self,gene_mutation=None,verbose=False,validate=True):
        '''
        Predict the effect of the given mutation on one or more antimicrobials.

        Args:
            mutation (str): a genetic variant in the form GENE_MUTATION e.g. for a SNP katG_S315T, INDEL katG_315_indel.
            verbose (bool): if True, then a description of the rules that apply to the supplied mutation and their priority is written to STDOUT (default=False)
            validate (bool): if True, then the supplied mutation is validated against the GenBank file. Only set to False if the mutation is being supplied from a trusted source and you know what you are doing! (default=True)

        Returns:
            result (dict): the drugs affected by the mutation are the keys, and the predicted phenotypes are the values. e.g. {'LEV':'R', 'MXF':'R'}
                           if the gene isn't in the catalogue, then an "S" is returned, on the assumption that it is probably susceptible.

        Notes:
            * mutations can be specified in a grammar that covers most of the known and expected genetic variants.
            * Stop codon is represented by "!"
            * "any mutation at position S315" (i.e. a wildcard) is represented by "?" e.g. S315?
            * for more info see the walkthrough and also the NOMENCLATURE.md file
        '''

        # split the supplied mutation on _ which separates the components
        components=gene_mutation.split("_")

        # the gene name is always the first component
        gene_name=components[0]

        # ..and the remainder is the mutation
        mutation=gene_mutation.split(gene_name+"_")[1]


        if validate:
            # check that the gene exists in the reference genome!
            assert self.reference_genome.valid_gene(gene_name), gene_name+" does not exist in the specified GENBANK file!"

            # also check that the supplied mutation is valid given the reference genome
            assert self.reference_genome.valid_mutation(gene_mutation), "gene exists but "+gene_mutation+" is badly formed; check the reference amino acid or nucleotide!"

        # parse the mutation to work out what type of mutation it is, where it acts etc
        (position,variant_affects,variant_type,indel_type,indel_length,indel_bases,before,after)=self._parse_mutation(gene_mutation)

        # insist we've been given an amino acid or a wildcard only
        if variant_type=="SNP":
            assert after in ['a','c','t','g','x','z','?',"!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid or base!"
            assert before in ['a','c','t','g','x','z','?','!','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], before+" is not an amino acid or base!"

        # check the position isn't lower than -100
        if variant_affects=="INDEL" and position!="*":
            assert position>=-100, "promoter mutation too far from CDS; max allowed is -100, position is "+str(position)

        if gene_name not in self.gene_list:

            if verbose:
                print(gene_name + " not in catalogue so assuming mutation is susceptible")

            # if it isn't, then we don't know what drugs it might affect, so just return susceptible
            return("S")

        # otherwise it _must_ be in our catalogue
        else:

            # find out the drugs to which changes in this gene can confer resistance
            drugs=self.gene_lookup[gene_name]

            # if there aren't any just return an S again
            if drugs is None:

                if verbose:
                    print(gene_name+" has no drugs associated with it in the catalogue, so assuming susceptible")

                return("S")

            # create the result dictionary e.g. {"RIF":'R', "RFB":'R'}
            # most of the time the predictions will be identical, but this allows them to diverge in future
            result={}

            position_vector=(self.resistance_catalogue.POSITION.isin([str(position),'*']))
            variant_affects_vector=(self.resistance_catalogue.VARIANT_AFFECTS==variant_affects)
            variant_type_vector=(self.resistance_catalogue.VARIANT_TYPE==variant_type)
            gene_vector=(self.resistance_catalogue.GENE==gene_name)
            # prediction_vector=(self.resistance_catalogue[self.catalogue_name+"_PREDICTION"].isin(["R","S","U"]))

            # deal with each compound, one at a time
            for compound in drugs:

                if verbose:
                    print(compound, gene_mutation)

                rules=self.resistance_catalogue.loc[ gene_vector &
                                                     position_vector &
                                                    (self.resistance_catalogue.DRUG==compound) &
                                                    variant_affects_vector &
                                                    variant_type_vector]

                # prepare a dictionary to store hits with the priority as the key:  e.g. {10:'R',5:'U'}
                predictions={}

                if not rules.empty:

                    rules_variant_type_vector=(rules.VARIANT_TYPE==variant_type)
                    rules_position_vector=(rules.POSITION==str(position))

                    if variant_type=="SNP":

                        # PRIORITY=1: synonymous mutation at any position in the CDS (e.g. rpoB_*=)
                        row=rules.loc[rules_variant_type_vector & (rules.MUTATION=="*=") & (before==after)]
                        if not row.empty:
                            predictions[1]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                            if verbose:
                                print("1. "+predictions[1]+" syn SNP at any position in the CDS")

                        # PRIORITY=10: nonsynoymous mutation at any position in the CDS or PROM (e.g. rpoB_*? or rpoB_-*?)
                        if variant_affects in ["CDS","RNA"]:
                            row=rules.loc[rules_variant_type_vector & (rules.MUTATION=="*?") & (before!=after)]
                        else:
                            row=rules.loc[rules_variant_type_vector & (rules.MUTATION=="-*?") & (before!=after)]
                        if not row.empty:
                            predictions[2]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                            if verbose:
                                print("2. "+predictions[2]+" nonsyn SNP at any position in the CDS or PROM")

                        # PRIORITY=100: any nonsynoymous mutation at this specific position in the CDS or PROM  (e.g. rpoB_S450? or rpoB_c-15?)
                        row=rules.loc[rules_variant_type_vector & rules_position_vector & (before!=after) & (rules.MUTATION.str[-1]=="?")]
                        if not row.empty:
                            predictions[3]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                            if verbose:
                                print("3. "+predictions[3]+" nonsyn SNP at specified position in the CDS")

                        # PRIORITY=1000: an exact match
                        row=rules.loc[rules_variant_type_vector & (rules.MUTATION==mutation)]
                        if not row.empty:
                            predictions[4]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                            if verbose:
                                print("4. "+predictions[4]+" exact SNP match")

                    elif variant_type=="INDEL":

                        # PRIORITY 1: any insertion or deletion in the CDS or PROM (e.g. rpoB_*_indel or rpoB_-*_indel)
                        if variant_affects in ["CDS","RNA"]:
                            row=rules.loc[rules_variant_type_vector & (rules.MUTATION=="*_indel")]
                        else:
                            row=rules.loc[rules_variant_type_vector & (rules.MUTATION=="-*_indel")]
                        if not row.empty:
                            predictions[1]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                            if verbose:
                                print("1. "+predictions[1]+" any insertion or deletion in the CDS or PROM")

                        # PRIORITY 2: rpoB_*_ins, rpoB_*_del any insertion (or deletion) in the CDS or PROM
                        if variant_affects in ["CDS","RNA"]:
                            row=rules.loc[rules_variant_type_vector & (rules.MUTATION.isin(["*_ins","*_del"]))]
                        else:
                            row=rules.loc[rules_variant_type_vector & (rules.MUTATION.isin(["-*_ins","-*_del"]))]
                        if not row.empty:
                            predictions[2]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                            if verbose:
                                print("2. "+predictions[2]+" any insertion (or deletion) in the CDS or PROM")

                        # PRIORITY 3: any insertion of a specified length (or deletion) in the CDS or PROM (e.g. rpoB_*_ins_2, rpoB_*_del_3, rpoB_-*_ins_1, rpoB_-*_del_200)
                        if indel_length is not None and indel_type!="indel":
                            if variant_affects in ["CDS","RNA"]:
                                row=rules.loc[rules_variant_type_vector & (rules.MUTATION.isin(["*_"+indel_type+"_"+str(indel_length)]))]
                            else:
                                row=rules.loc[rules_variant_type_vector & (rules.MUTATION.isin(["-*_"+indel_type+"_"+str(indel_length)]))]
                            if not row.empty:
                                predictions[3]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                                if verbose:
                                    print("3. "+predictions[3]+" any insertion of a specified length (or deletion) in the CDS or PROM")

                        # PRIORITY=4: any frameshifting insertion or deletion in the CDS (e.g. rpoB_*_fs)
                        if indel_length is not None and (indel_length%3)==0:
                            row=rules.loc[rules_variant_type_vector & (rules.MUTATION=="*_fs")]
                            if not row.empty:
                                predictions[4]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                                if verbose:
                                    print("4. "+predictions[4]+" any frameshifting insertion or deletion in the CDS")

                        # PRIORITY=5: any indel at a specific position in the CDS or PROM (e.g. rpoB_1300_indel or rpoB_-15_indel)
                        row=rules.loc[rules_variant_type_vector & rules_position_vector & (rules.MUTATION.str.contains("indel"))]
                        if not row.empty:
                            predictions[5]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                            if verbose:
                                print("5. "+predictions[5]+" any indel at a specific position in the CDS or PROM (e.g. rpoB_1300_indel or rpoB_-15_indel)")

                        # PRIORITY=6: an insertion (or deletion) at a specific position in the CDS or PROM (e.g. rpoB_1300_ins or rpoB_-15_del)
                        if indel_type!="indel":
                            row=rules.loc[rules_variant_type_vector & rules_position_vector & (rules.MUTATION.str.contains(indel_type))]
                            if not row.empty:
                                predictions[6]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                                if verbose:
                                    print("6. "+predictions[6]+" any insertion (or deletion) at a specific position in the CDS or PROM (e.g. rpoB_1300_ins or rpoB_-15_del)")

                        # PRIORITY=7: an insertion (or deletion) of a specified length at a specific position in the CDS or PROM (e.g. rpoB_1300_ins_2 or rpoB_-15_del_200)
                        if indel_type!="indel" and indel_length is not None:
                            row=rules.loc[rules_variant_type_vector & rules_position_vector & (rules.MUTATION==str(position)+"_"+indel_type+"_"+str(indel_length))]
                            if not row.empty:
                                predictions[7]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                                if verbose:
                                    print("7. "+predictions[7]+" an insertion (or deletion) of a specified length at a specific position in the CDS or PROM (e.g. rpoB_1300_ins_2 or rpoB_-15_del_200)")

                        # PRIORITY=8: a frameshifting mutation at a specific position in the CDS (e.g. rpoB_100_fs)
                        if indel_length is not None and (indel_length%3)==0 and position is not None:
                            row=rules.loc[rules_variant_type_vector & rules_position_vector & (rules.MUTATION==str(position)+"_fs")]
                            if not row.empty:
                                predictions[8]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                                if verbose:
                                    print("8. "+predictions[8]+" a frameshifting mutation at a specific position in the CDS (e.g. rpoB_100_fs)")

                        # PRIORITY=9: an insertion of a specified series of nucleotides at a position in the CDS or PROM (e.g. rpoB_1300_ins_ca)
                        if indel_type!="indel" and indel_length is not None and indel_bases is not None:
                            row=rules.loc[rules_variant_type_vector & rules_position_vector & (rules.MUTATION==str(position)+"_"+indel_type+"_"+str(indel_length)+"_"+indel_bases)]
                            if not row.empty:
                                predictions[9]=str(row[self.catalogue_name+"_PREDICTION"].values[0])
                                if verbose:
                                    print("9. "+predictions[9]+" an insertion of a specified series of nucleotides at a position in the CDS or PROM (e.g. rpoB_1300_ins_ca)")
                    if verbose:
                        print("----------------------")


                if not predictions:
                    # all mutations should hit at least one of the default entries, so if this doesn't happen, something is wrong
                    # return({"UNK":"U"})
                    raise ValueError("No entry found in the catalogue for "+gene_mutation)
                else:
                    # if the dictionary is not empty, store the result with the highest priority
                    final_prediction=predictions[sorted(predictions)[-1]]

                result[compound]=final_prediction

        return(result)

    def _parse_mutation(self,gene_mutation):
        """
        Take a GENE_MUTATION, and determine what type it is, what it affects etc

        Args:
            gene_mutation (str): the form as defined by GOARC e.g. rpoB_S450L or rpoB_1300_ins_3

        Returns:
            position (int): the position of the mutation. This is the amino acid number in a protein, or the number of the nucleotide in the promoter
            variant_affects (str): whether this genetic variant affects the promoter (PROM), coding sequence (CDS)
            variant_type (str): SNP or INDEL
            indel_type (str): indel, ins, del or fs
            indel_length (int): if sufficient information is given, the number of bases in the INDEL
            indel_bases (str): if sufficient information is given, the bases in the INDEL
            before (str): the REF amino acid
            after (str): the ALT amino acid

        """

        variant_type=None
        variant_affects=None
        position=None
        indel_type=None
        indel_length=None
        indel_bases=None
        before=None
        after=None

        cols=gene_mutation.split("_")

        if len(cols)==2:

            variant_type="SNP"

            mutation=cols[1]

            position=int(mutation[1:-1])
            if position<0:
                variant_affects="PROM"
            else:
                if self.reference_genome.gene_type(cols[0])=='RNA':
                    variant_affects="RNA"
                else:
                    variant_affects="CDS"

            before=mutation[0]
            after=mutation[-1]

        elif len(cols) in [3,4]:

            variant_type="INDEL"
            before=None
            after=None

            position=int(cols[1])
            if position<0:
                variant_affects="PROM"
            else:
                if self.reference_genome.gene_type(cols[0])=='RNA':
                    variant_affects="RNA"
                else:
                    variant_affects="CDS"

            # the third element is one of indel, ins, del or the special case fs
            indel_type=cols[2]

            # if there is a fourth and final element to an INDEL it is either _4 or _ctgc
            if len(cols)==4:
                if cols[3].isnumeric():
                    indel_length=int(cols[3])
                    indel_bases=None
                else:
                    indel_length=len(cols[3])
                    indel_bases=cols[3]

        return(position,variant_affects,variant_type,indel_type,indel_length,indel_bases,before,after)
