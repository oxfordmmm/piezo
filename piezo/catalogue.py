#! /usr/bin/env python

import os
import collections

import pandas, numpy, ujson


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
                                      'rules' ]
                                    )

class ResistanceCatalogue:

    def __init__(self, catalogue_file):

        self.catalogue = load_catalogue(catalogue_file)

    def predict(self, mutation, verbose=False):

        return predict(self.catalogue, mutation=mutation, verbose=verbose)

def parse_json(data):

    return ujson.loads(data)

def split_mutation(row):

    # split the supplied mutation on _ which separates the components
    components=row['MUTATION'].split("_")

    # the gene name is always the first component
    gene=components[0]

    # ..and the remainder is the mutation
    mutation=row['MUTATION'].split(gene+"_")[1]

    if len(components)==2:
        mutation_type="SNP"
        mutation_affects="CDS"
        try:
            position=components[1][1:-1]
            if int(position)<0:
                mutation_affects="PROM"
            elif int(position)==0:
                raise ValueError("a position of 0 in a mutation does not make sense!")
        except:
            position=components[1][:-1]
            if position[0]=="-":
                mutation_affects="PROM"

    elif len(components) in [3,4]:
        mutation_type="INDEL"
        mutation_affects="CDS"
        position=components[1]
        try:
            if int(position)<0:
                mutation_affects="PROM"
            elif int(position)==0:
                raise ValueError("a position of 0 in a mutation does not make sense!")
        except:
            try:
                if position[0]=="-":
                    mutation_affects="PROM"
            except:
                print(position)

    else:
        raise ValueError("badly formed mutation: "+row["MUTATION"])

    return(pandas.Series([gene,mutation,position,mutation_affects,mutation_type]))


def load_catalogue(catalogue_file):

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

    assert len(rules.CATALOGUE_VALUES.unique())==1, "multiple grammar value types used in the catalogue!"
    values=rules.CATALOGUE_VALUES.unique()[0]
    if values in ['RS','RUS','RFUS']:
        values=[i for i in rules.CATALOGUE_VALUES.unique()[0]]
    else:
        raise ValueError("content of column CATALOGUE_VALUES not recognised!")

    # since we are storing the above values in the named tuple, we don't need them in the rules dataframe
    rules=rules[['DRUG','MUTATION','PREDICTION','EVIDENCE','SOURCE','OTHER']]

    # find out the DRUGS that this catalogue applies to
    drugs=list(rules.DRUG.unique())

    if grammar!="GM1":

        raise ValueError("only the GENE_MUTATION grammar is supported at present!")

    elif grammar=="GM1":

        rules[['GENE','MUTATION','POSITION','MUTATION_AFFECTS','MUTATION_TYPE']]=rules.apply(split_mutation,axis=1)

        genes=list(rules.GENE.unique())

        gene_lookup={}
        drug_lookup={}

        for gene_name in genes:

            df=rules.loc[rules.GENE==gene_name]

            gene_lookup[gene_name]=list(df.DRUG.unique())

        for drug_name in drugs:

            df=rules.loc[rules.DRUG==drug_name]

            drug_lookup[drug_name]=list(df.GENE.unique())


        return catalogue(genbank_reference,
                               name,
                               version,
                               grammar,
                               values,
                               drugs,
                               genes,
                               drug_lookup,
                               gene_lookup,
                               rules)

def predict(catalogue,mutation,verbose):

    if catalogue.grammar=="GM1":

        result = predict_GM1(catalogue,mutation,verbose)

        return(result)

    else:

        raise ValueError("Only the GENE_MUTATION grammar is supported at present")


def predict_GM1(catalogue,gene_mutation,verbose):

    components=gene_mutation.split("_")

    gene=components[0]
    mutation=gene_mutation.split(gene+"_")[1]

    # parse the mutation to work out what type of mutation it is, where it acts etc
    (position, mutation_affects, mutation_type, indel_type, indel_length, indel_bases, before, after) = parse_GM1(gene,mutation)

    # insist we've been given an amino acid or a wildcard only
    if mutation_type=="SNP":
        assert after in ['a','c','t','g','x','z','?',"!",'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not an amino acid or base!"
        assert before in ['a','c','t','g','x','z','?','!','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z'], before+" is not an amino acid or base!"

    # if the gene isn't in the catalogue, then assume it has no effect
    if gene not in catalogue.genes:
        return("S")

    # find out the drugs to which changes in this gene can confer resistance
    drugs=catalogue.gene_lookup[gene]

    # if there aren't any just return an S again
    if drugs is None:
        return("S")

    # create the result dictionary e.g. {"RIF":'R', "RFB":'R'}
    # most of the time the predictions will be identical, but this allows them to diverge in future
    result={}

    # create the vectors of Booleans that will apply
    position_vector = (catalogue.rules.POSITION.isin([position,str(position), '*']))
    mutation_affects_vector = (catalogue.rules.MUTATION_AFFECTS == mutation_affects)
    mutation_type_vector = (catalogue.rules.MUTATION_TYPE == mutation_type)
    gene_vector = (catalogue.rules.GENE == gene)

    rules_subset_vector = position_vector & mutation_affects_vector & mutation_type_vector & gene_vector

    # deal with each compound, one at a time
    for compound in drugs:

        rules_subset = catalogue.rules.loc[(rules_subset_vector) & (catalogue.rules.DRUG==compound)]

        # prepare a dictionary to store hits with the priority as the key:  e.g. {10:'R',5:'U'}
        predictions={}

        if not rules_subset.empty:

            position_subset_vector=(rules_subset.POSITION==str(position))
            mutation_type_subset_vector=(rules_subset.MUTATION_TYPE==mutation_type)

            if mutation_type=="SNP":
                process_snp_variants(mutation_affects, predictions, catalogue, before, after, mutation, rules_subset,
                                     mutation_type_subset_vector, position_subset_vector, verbose)
            elif mutation_type=="INDEL":
                process_indel_variants(variant_affects,
                       predictions, catalogue, before, after, mutation, rules_subset, rules_mutation_type_vector,
                       rules_position_vector, indel_length, indel_type, indel_bases, position, verbose)

        if not predictions:
            # all mutations should hit at least one of the default entries, so if this doesn't happen, something is wrong
            raise ValueError("No entry found in the catalogue for "+gene_mutation)

        final_prediction=predictions[sorted(predictions)[-1]]
        result[compound]=final_prediction

    return(result)


def row_prediction(row, predictions, priority, message):
    if not row.empty:
        predictions[priority] = str(row["PREDICTION"].values[0])


def process_snp_variants(mutation_affects,
                         predictions,
                         catalogue,
                         before, after,
                         mutation,
                         rules,
                         rules_mutation_type_vector,
                         rules_position_vector,
                         verbose):
    # PRIORITY=1: synonymous mutation at any position in the CDS (e.g. rpoB_*=)
    row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*=") & (before==after)]
    row_prediction(row, predictions, 1, "syn SNP at any position in the CDS")

    # PRIORITY=2: nonsynoymous mutation at any position in the CDS or PROM (e.g. rpoB_*? or rpoB_-*?)
    if mutation_affects=="CDS":
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*?") & (before!=after)]
    else:
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="-*?") & (before!=after)]
    row_prediction(row, predictions, 2, "nonsyn SNP at any position in the CDS or PROM")

    # PRIORITY=3: het mutation at any position in the CDS or PROM (e.g. rpoB_*Z or rpoB_-*z)
    if mutation_affects=="CDS":
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*Z") & (before!=after) & (mutation[-1]=="Z")]
    else:
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="-*z") & (before!=after) & (mutation[-1]=="z")]
    row_prediction(row, predictions, 3, "het SNP at any position in the CDS or PROM")

    # PRIORITY=4: any nonsynoymous mutation at this specific position in the CDS or PROM  (e.g. rpoB_S450? or rpoB_c-15?)
    row=rules.loc[rules_mutation_type_vector & rules_position_vector & (before!=after) & (rules.MUTATION.str[-1]=="?")]
    row_prediction(row, predictions, 4, "nonsyn SNP at specified position in the CDS")

    # PRIORITY=5: any het mutation at this specific position in the CDS or PROM  (e.g. rpoB_S450Z or rpoB_c-15z)
    row=rules.loc[(rules_mutation_type_vector) & (rules_position_vector) & (before!=after) & (mutation[-1] in ["Z","z"]) & (rules.MUTATION.str[-1] in ['Z','z'])]
    row_prediction(row, predictions, 5, "het SNP at specified position in the CDS")

    # PRIORITY=6: an exact match
    row=rules.loc[(rules_mutation_type_vector) & (rules.MUTATION==mutation)]
    row_prediction(row, predictions, 6, "exact SNP match")

def process_indel_variants(mutation_affects,
                           predictions,
                           catalogue,
                           before, after,
                           mutation,
                           rules,
                           rules_mutation_type_vector,
                           rules_position_vector,
                           indel_length,
                           indel_type,
                           indel_bases,
                           position,
                           verbose):
    # PRIORITY 1: any insertion or deletion in the CDS or PROM (e.g. rpoB_*_indel or rpoB_-*_indel)
    if mutation_affects == "CDS":
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*_indel")]
    else:
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="-*_indel")]
    row_prediction(row, predictions, 1, "any insertion or deletion in the CDS or PROM")

    # PRIORITY 2: rpoB_*_ins, rpoB_*_del any insertion (or deletion) in the CDS or PROM
    if mutation_affects == "CDS":
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION.isin(["*_ins","*_del"]))]
    else:
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION.isin(["-*_ins","-*_del"]))]
    row_prediction(row, predictions, 2, "any insertion (or deletion) in the CDS or PROM")

    # PRIORITY 3: any insertion of a specified length (or deletion) in the CDS or PROM (e.g. rpoB_*_ins_2, rpoB_*_del_3, rpoB_-*_ins_1, rpoB_-*_del_200)
    if indel_length is not None and indel_type!="indel":
        if mutation_affects == "CDS":
            row=rules.loc[rules_mutation_type_vector & (rules.MUTATION.isin(["*_"+indel_type+"_"+str(indel_length)]))]
        else:
            row=rules.loc[rules_mutation_type_vector & (rules.MUTATION.isin(["-*_"+indel_type+"_"+str(indel_length)]))]
        row_prediction(row, predictions, 3, "any insertion of a specified length (or deletion) in the CDS or PROM")

    # PRIORITY=4: any frameshifting insertion or deletion in the CDS (e.g. rpoB_*_fs)
    if indel_length is not None and (indel_length%3)==0:
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*_fs")]
        row_prediction(row, predictions, 4, "any frameshifting insertion or deletion in the CDS")

    # PRIORITY=5: any indel at a specific position in the CDS or PROM (e.g. rpoB_1300_indel or rpoB_-15_indel)
    row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION.str.contains("indel"))]
    row_prediction(row, predictions, 5, "any indel at a specific position in the CDS or PROM (e.g. rpoB_1300_indel or rpoB_-15_indel)")

    # PRIORITY=6: an insertion (or deletion) at a specific position in the CDS or PROM (e.g. rpoB_1300_ins or rpoB_-15_del)
    if indel_type!="indel":
        row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION.str.contains(indel_type))]
        row_prediction(row, predictions, 6, "any insertion (or deletion) at a specific position in the CDS or PROM (e.g. rpoB_1300_ins or rpoB_-15_del)")

    # PRIORITY=7: an insertion (or deletion) of a specified length at a specific position in the CDS or PROM (e.g. rpoB_1300_ins_2 or rpoB_-15_del_200)
    if indel_type!="indel" and indel_length is not None:
        row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION==str(position)+"_"+indel_type+"_"+str(indel_length))]
        row_prediction(row, predictions, 7, "an insertion (or deletion) of a specified length at a specific position in the CDS or PROM (e.g. rpoB_1300_ins_2 or rpoB_-15_del_200)")

    # PRIORITY=8: a frameshifting mutation at a specific position in the CDS (e.g. rpoB_100_fs)
    if indel_length is not None and (indel_length%3)==0 and position is not None:
        row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION==str(position)+"_fs")]
        row_prediction(row, predictions, 8, "a frameshifting mutation at a specific position in the CDS (e.g. rpoB_100_fs)")

    # PRIORITY=9: an insertion of a specified series of nucleotides at a position in the CDS or PROM (e.g. rpoB_1300_ins_ca)
    if indel_type!="indel" and indel_length is not None and indel_bases is not None:
        row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION==str(position)+"_ins_"+indel_bases)]
        row_prediction(row, predictions, 9, "an insertion of a specified series of nucleotides at a position in the CDS or PROM (e.g. rpoB_1300_ins_ca)")


def parse_GM1(gene,mutation):

    mutation_type=None
    mutation_affects=None
    position=None
    indel_type=None
    indel_length=None
    indel_bases=None
    before=None
    after=None

    cols=mutation.split("_")

    if len(cols)==1:

        mutation_type="SNP"

        position=int(mutation[1:-1])
        if position<0:
            mutation_affects="PROM"
        else:
            mutation_affects="CDS"

        before=mutation[0]
        after=mutation[-1]

    elif len(cols) in [2,3]:
        mutation_type="INDEL"

        position=int(cols[0])
        if position<0:
            mutation_affects="PROM"
        else:
            mutation_affects="CDS"

        # the third element is one of indel, ins, del or the special case fs
        indel_type=cols[1]

        # if there is a fourth and final element to an INDEL it is either _4 or _ctgc
        if len(cols)==3:
            if cols[2].isnumeric():
                indel_length=int(cols[2])
            else:
                indel_length=len(cols[2])
                indel_bases=cols[3]

    return(position, mutation_affects, mutation_type, indel_type, indel_length, indel_bases, before, after)
