#! /usr/bin/env python
'''All logic required to parse and generate predictions from a catalogue in GARC1
'''

import collections
import re

import pandas


def split_mutation(row: pandas.Series) -> pandas.Series:
    '''Take a row of the catalogue and get info about the mutation from it.
    This includes the type, ref, alt, position and minor populations

    Args:
        row (pandas.Series): Catalogue row

    Raises:
        ValueError: Raised when mutations are malformed

    Returns:
        pandas.Series: Series about the mutation. In order: gene, mutation, position, type of position (i.e PROM/CDS), mutation type (i.e INDEL/SNP/MULTI), support for minor populations
    '''
    if '&' in row['MUTATION']:
        #Multi-mutation so don't try to split it out in the same way
        gene = 'MULTI'
        #Ensure mutations are in a reproducable order
        mutation = '&'.join(sorted(list(row['MUTATION'].split('&'))))

        position = None
        mutation_affects = None
        mutation_type = 'MULTI'

        #Check for minor populations (and remove, we'll check these separately)
        minors = ','.join([mut.split(":")[-1] if ':' in mut else '' for mut in mutation.split("&")])
        mutation = '&'.join([mut.split(":")[0] for mut in mutation.split("&")])

        return pandas.Series([gene, mutation, position, mutation_affects, mutation_type, minors])

    # split the supplied mutation on _ which separates the components
    components=row['MUTATION'].split("@")

    # the gene name is always the first component
    gene=components[0]

    # ..and the remainder is the mutation
    mutation=row['MUTATION'].split(gene+"@")[1]

    #Check for minor populations (and remove minors, we'll check for those separately)
    if ":" in mutation:
        minor = mutation.split(":")[-1]
        mutation = mutation.split(":")[0]
    else:
        minor = ''


    cols=mutation.split('_')
    valid = True

    if len(cols)==1:
        mutation_type="SNP"
        mutation_affects="CDS"
        try:
            position=cols[0][1:-1]
            if int(position)<0:
                mutation_affects="PROM"
            elif int(position)==0:
                valid = False
        except:
            position=cols[0][:-1]
            if position[0]=="-":
                mutation_affects="PROM"

    elif len(cols) in [2,3]:
        mutation_type="INDEL"
        percentage = False
        if len(cols) == 2:
            #Checking for percentage deletions (of form "del_<percent>")
            try:
                float(cols[1])
                percentage = True
            except ValueError:
                percentage = False
            if percentage:
                position = None
        if percentage:
            mutation_affects = "GENE"
        else:
            mutation_affects = "CDS"
            position=cols[0]
            try:
                if int(position)<0:
                    mutation_affects="PROM"
                elif int(position)==0:
                    valid = False
            except:
                try:
                    if position[0]=="-":
                        mutation_affects="PROM"
                except:
                    valid = False

    else:
        raise ValueError("Badly formed mutation: "+row["MUTATION"])
    
    if not valid:
        #Valid mutation not found, so complain
        raise ValueError("Badly formed mutation: "+row['MUTATION'])
    
    return pandas.Series([gene,mutation,position,mutation_affects,mutation_type, minor])


def process_catalogue_GARC1(rules: pandas.DataFrame, drugs: [str], catalogue_genes_only: bool) -> (pandas.DataFrame, [str], {str: [str]}, {str: [str]}):
    '''
    For the GARC1 grammar, add some additional columns to the rules dataframe that can be inferred from the mutation.

    Args:
        rules (pandas.DataFrame): Pandas dataframe of the AMR catalogue in the general form
        drugs ([str]): list of drugs in the AMR catalogue
        catalogue_genes_only (bool): whether to subset the catalogue down so it ONLY includes (DRUG,GENE) pairs that include at least one row predicting resistance
    Returns:
        (pandas.DataFrame, [str], {str: [str]}, {str: [str]}) : Tuple of values. In order:
                                                                    DataFrame of the rules
                                                                    List of the genes
                                                                    Dictionary mapping drug name -> [gene names]
                                                                    Dictionary mapping gene name -> [drug names]

    '''

    # infer these columns
    rules[['GENE','MUTATION','POSITION','MUTATION_AFFECTS','MUTATION_TYPE', "MINOR"]] = rules.apply(split_mutation,axis=1)

    if catalogue_genes_only:

        # create an index of which (drug,gene) pairs have at least one row predicting resistance or uncertainty
        drug_list_index=rules.loc[(rules.PREDICTION=="R") | (rules.PREDICTION == 'U')][['DRUG','GENE']].groupby(['DRUG','GENE']).count()

        # index the rules so we can right join
        rules.set_index(['DRUG','GENE'],inplace=True)

        # do the right join; this has the effect of removing all the rows for "genes of interest"
        rules=rules.join(drug_list_index,how='right')

        # remove the index
        rules.reset_index(inplace=True)

    # create a list of the genes mentioned in the catalogue
    genes=list(rules.GENE.unique())

    # create a dictionary where the genes are keys and the entries tell which drugs are affected
    gene_lookup={}
    for gene_name in genes:
        df=rules.loc[rules.GENE==gene_name]
        gene_lookup[gene_name] = sorted(list(df.DRUG.unique()))

    # create a dictionary where the drugs are keys and the entries tell which genes are affected
    drug_lookup={}
    for drug_name in drugs:
        df=rules.loc[rules.DRUG==drug_name]
        drug_lookup[drug_name]=list(df.GENE.unique())

    return rules, genes, drug_lookup, gene_lookup

def predict_GARC1(catalogue: collections.namedtuple, gene_mutation: str) -> {str: str}:
    '''
    For the GARC1 grammar, predict the effect of the supplied gene_mutation.

    Args:
        catalogue (collections.namedtuple): defines the resistance catalogue
        gene_mutation (str): as specified by the GARC1 grammar, e.g. rpoB_S450L, fabG1_c-15t, ahpC_128_indel

    Returns:
        {str: str} | str: Either a dictionary mapping drug name -> predictions or "S" if no hits in the catalogue
    '''
    if '&' in gene_mutation:
        #Multi-gene so treat differently
        return predict_multi(catalogue, gene_mutation)
    components=gene_mutation.split("@")

    gene=components[0]
    mutation=gene_mutation.split(gene+"@")[1]

    #Check for minors
    if ":" in mutation:
        try:
            minor = float(mutation.split(":")[-1])
        except ValueError:
            assert False, "Malformed mutation! "+mutation
        mutation = mutation.split(":")[0]
        assert minor > 0, "Minor population given with no evidence! "+mutation
    else:
        minor = None

    # parse the mutation to work out what type of mutation it is, where it acts etc
    (position, mutation_affects, mutation_type, indel_type, indel_length, indel_bases, before, after) = parse_mutation(mutation)

    # if the gene isn't in the catalogue, then assume it has no effect
    if gene not in catalogue.genes:
        return("S")

    # find out the drugs to which changes in this gene can confer resistance
    drugs=catalogue.gene_lookup[gene]

    # create the result dictionary e.g. {"RIF":'R', "RFB":'R'}
    # most of the time the predictions will be identical, but this allows them to diverge in future
    result={}


    # create the vectors of Booleans that will apply
    position_vector = (catalogue.rules.POSITION.isin([position,str(position), '*','-*']))
    mutation_affects_vector = (catalogue.rules.MUTATION_AFFECTS == mutation_affects)
    mutation_type_vector = (catalogue.rules.MUTATION_TYPE == mutation_type)
    gene_vector = (catalogue.rules.GENE == gene)
    subset_vector = position_vector & mutation_affects_vector & mutation_type_vector & gene_vector
    # subset_vector = (catalogue.rules.POSITION.isin([position,str(position), '*','-*'])) & (catalogue.rules.MUTATION_AFFECTS == mutation_affects) & (catalogue.rules.MUTATION_TYPE == mutation_type) & (catalogue.rules.GENE == gene)
    # deal with each compound, one at a time
    for compound in drugs:

        subset_rules = catalogue.rules.loc[(subset_vector) & (catalogue.rules.DRUG==compound)]

        # prepare a dictionary to store hits with the priority as the key:  e.g. {10:'R',5:'U'}
        predictions={}

        if not subset_rules.empty:

            subset_position_vector=(subset_rules.POSITION==str(position))
            subset_mutation_type_vector=(subset_rules.MUTATION_TYPE==mutation_type)

            if mutation_type=="SNP":
                process_snp_variants(mutation_affects, predictions, before, after, mutation, subset_rules,\
                                    subset_mutation_type_vector, subset_position_vector, minor)
            elif mutation_type=="INDEL":
                process_indel_variants(mutation_affects,
                       predictions, subset_rules, subset_mutation_type_vector,
                       subset_position_vector, indel_length, indel_type, indel_bases, position, minor)

        if not predictions:
            # all mutations should hit at least one of the default entries, so if this doesn't happen, something is wrong
            raise ValueError("No entry found in the catalogue for "+gene_mutation+" "+compound)

        final_prediction=predictions[sorted(predictions)[-1]]
        result[compound]=final_prediction

    return(result)

def predict_multi(catalogue: collections.namedtuple, gene_mutation: str) -> {str: str}:
    '''Get the predictions for a given multi-mutation. 
    Mutli-mutations are (currently) a lot stricter than others, and do not support wildcards or subsets of the mutations

    Args:
        catalogue (collections.namedtuple): The resistance catalogue
        gene_mutation (str): String of the mutations, each in GARC, separated by `&`
    Returns:
        {str: str} | str : Either a dict with drugs as the keys, or 'S' when no predictions in catalogue
    '''
    #Ensure that the mutations are in a reproducable order
    sorted_mutation = '&'.join(sorted(list(gene_mutation.split('&'))))

    #Check for minor populations (and remove, we'll check these separately)
    minors = [float(mut.split(":")[-1]) if ':' in mut else None for mut in sorted_mutation.split("&")]
    sorted_mutation = '&'.join([mut.split(":")[0] for mut in sorted_mutation.split("&")])

    #Get the multi rules
    multi_rules = catalogue.rules[catalogue.rules['MUTATION_TYPE']=='MULTI']

    drugs = {drug: 'S' for drug in catalogue.drugs}
    values = catalogue.values
    for (_, rule) in multi_rules.iterrows():
        if rule['MUTATION'] == sorted_mutation:
            #We have a match! Prioritise predictions based on values
            drug = rule['DRUG']
            if values.index(rule['PREDICTION']) < values.index(drugs[drug]):
                #The prediction is closer to the start of the values list, so should take priority
                #Check for minor populations first though
                for (cat, minor) in zip(rule['MINOR'].split(","), minors):
                    if minor is None and cat == '':
                        #Neither are minors so add
                        drugs[drug] = rule['PREDICTION']
                    elif cat != '' and minor is not None and minor < 1 and float(cat) < 1:
                        #FRS
                        if minor >= float(cat):
                            #Match
                            drugs[drug] = rule['PREDICTION']
                    elif cat != '' and minor is not None and minor >=1 and float(cat) >= 1:
                        #COV
                        if minor >= float(cat):
                            #Match
                            drugs[drug] = rule['PREDICTION']

    #Check to ensure we have at least 1 prediction
    if len([key for key in drugs.keys() if drugs[key] != 'S']) > 0:
        return {drug: drugs[drug] for drug in drugs.keys() if drugs[drug]!='S'}

    #Nothing predicted, so try each individual mutation
    predictions = {}
    for (mutation, minor) in zip(sorted_mutation.split("&"), minors):
        #Put minor populations back into mutations as required
        if minor is not None:
            mutation = mutation + ":" + str(minor)

        #Get the prediction for the mutation
        pred = predict_GARC1(catalogue, mutation)
        if isinstance(pred, str):
            #We have a 'S' prediction so ignore it
            continue
        else:
            predictions[mutation] = pred

    return merge_predictions(predictions, catalogue)


def merge_predictions(predictions: {str: {str: str}}, catalogue: collections.namedtuple) -> {str: str}:
    '''When multi-mutations do not have a hit, they are decomposed and individuals are tried instead
    This merges these predictions based on the prioritisation defined for this catalogue

    Args:
        predictions ({str: {str: str}}): Dictionary mapping mutation -> effects for each mutation in the multi
        catalogue (collections.namedtuple): The catalogue object (for finding the values)
    Returns:
        {str: str}: Merged effects dict mapping drug name -> prediction. Or "S" if no matches
    '''
    #Pull out all of the drugs which have predictions
    drugs = set()
    for pred in predictions:
        for drug in predictions[pred].keys():
            drugs.add(drug)

    #Default to 'S' for now
    merged = {drug: 'S' for drug in drugs}

    #Look for all predictions for each drug, and report the most significant
    for drug in drugs:
        for mutation in predictions.keys():
            #Get the prediction for this drug (if exists)
            pred = predictions[mutation][drug]
            if catalogue.values.index(pred) <= catalogue.values.index(merged[drug]):
                #New highest priority, so set the prediction
                merged[drug] = pred
    
    #If we have >=1 non-'S' prediction, return the dict
    if len([key for key in merged.keys() if merged[key] != 'S']) > 0:
        return {drug: merged[drug] for drug in merged.keys() if merged[drug]!='S'}

    #Else give the default 'S'
    return 'S'

def row_prediction(rows: pandas.DataFrame, predictions: {int: str}, priority: int, minor: float) -> None:
    '''Get the predictions from the catalogue for the applicable rows

    Args:
        rows (pandas.DataFrame): DataFrame of the rows within the catalogue
        predictions (dict): Predictions dict mapping priorities to predictions. This is updated for implict return
        priority (int): Priority of this mutation type
        minor (float | None): Reads/FRS supporting this, or None if not a minor population
    '''
    pred = None
    values = ['R', 'U', 'F', 'S', None]
    for _, row in rows.iterrows():
        if not row.empty:
            assert int(priority) in range (1,11), 'priority must be an integer in range 1,2..10'
            if values.index(row['PREDICTION']) < values.index(pred):
                #This row's prediction is more important than the current, so check for minors and prioritise
                if minor is None and row['MINOR'] == '':
                    #Neither are minor populations so act normally
                    pred = row['PREDICTION']

                elif row['MINOR'] != '' and minor < 1 and float(row['MINOR']) < 1:
                    #We have FRS
                    if minor >= float(row['MINOR']):
                        #Match
                        pred = row['PREDICTION']

                elif row['MINOR'] != '' and minor >= 1 and float(row['MINOR']) >= 1:
                    #We have COV
                    if minor >= float(row['MINOR']):
                        #Match
                        pred = row['PREDICTION']

                else:
                    #Mismatch, so only give a match for defaults
                    if "*" in row['MUTATION']:
                        pred = row['PREDICTION']
    if pred:
        predictions[int(priority)] = str(pred)

def large_del(predictions: {int: str}, rules: pandas.DataFrame, size: float, minor: float) -> None:
    '''Row prediction, but specifically for large deletions

    Args:
        predictions ({int: str}): Predictions dictionary
        rules (pandas.DataFrame): Rules from the catalogue
        size (float): Percentage of the gene deleted (as a decimal)
        minor (float): Minors (if existing)
    '''
    pred = None
    current = -1
    #Find which rules are large deletions and see if the mutation is >= rule
    deletion = re.compile(r"""del_([01]\.[0-9]+)""")
    for _, rule in rules.iterrows():
        if deletion.fullmatch(rule['MUTATION']):
            percentage = float(deletion.fullmatch(rule['MUTATION']).groups()[0])
            if size >= percentage and percentage > current:
                #Match!
                if minor is None and rule['MINOR'] == '':
                    #Neither are minor populations so act normally
                    pred = rule['PREDICTION']
                    current = percentage
                elif minor is not None and rule['MINOR'] != '' and minor < 1 and float(rule['MINOR']) < 1:
                    #We have FRS
                    if minor >= float(rule['MINOR']):
                        #Match
                        pred = rule['PREDICTION']
                        current = percentage

                elif minor is not None and rule['MINOR'] != '' and minor >= 1 and float(rule['MINOR']) >= 1:
                    #We have COV
                    if minor >= float(rule['MINOR']):
                        #Match
                        pred = rule['PREDICTION']
                        current = percentage
                else:
                    #Mismatch so only match default
                    if rule['MUTATION'] == "del_0.0":
                        pred = rule['PREDICTION']
                        current = percentage
    if pred:
        predictions[1] = str(pred)


def process_snp_variants(mutation_affects: str,
                         predictions: {int: str},
                         before: str, after: str,
                         mutation: str,
                         rules: pandas.DataFrame,
                         rules_mutation_type_vector: pandas.Series,
                         rules_position_vector: pandas.Series,
                         minor: float) -> None:
    '''Apply cascading rules for SNPs according to GARC

    Args:
        mutation_affects (str): Where this mutation affects. i.e PROM/CDS
        predictions ({int: str}): Predictions dictionary mapping priority -> prediction. This is updated for implict return
        before (str): Reference base/AA
        after (str): Alt base/AA
        mutation (str): Mutation in GARC
        rules (pandas.DataFrame): Rules DataFrame
        rules_mutation_type_vector (pandas.Series): Series relating to rules for this type of mutation
        rules_position_vector (pandas.Series): Series relating to rules for this position
        minor (float | None): Float for supporting evidence of minor populations (or None if not a minor population).
    '''


    if before==after:

        # PRIORITY=0: no change, i.e. wildtype specified

        # PRIORITY=1: synonymous mutation at any position in the CDS (e.g. rpoB_*=)
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*=")]
        #syn SNP at any position in the CDS
        row_prediction(row, predictions, 1, minor)

    elif before!=after:

        # PRIORITY=2: nonsynoymous mutation at any position in the CDS or PROM (e.g. rpoB_*? or rpoB_-*?)
        if mutation_affects=="CDS":
            row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*?")]
            #nonsyn SNP at any position in the CDS
            row_prediction(row, predictions, 2, minor)
        else:
            row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="-*?")]
            #nonsyn SNP at any position in the PROM
            row_prediction(row, predictions, 2, minor)
        # PRIORITY=3: het mutation at any position in the CDS or PROM (e.g. rpoB_*Z or rpoB_-*z)
        if mutation[-1] in ['Z','z','O','o','X','x']:
            if mutation_affects=="CDS":
                if (mutation[-1]=="Z"):
                    row=rules.loc[(rules_mutation_type_vector) & (rules.MUTATION=="*Z")]
                    #het SNP at any position in the CDS
                    row_prediction(row, predictions, 3, minor)
                elif (mutation[-1]=="O"):
                    row=rules.loc[(rules_mutation_type_vector) & (rules.MUTATION=="*O")]
                    #filter fail at any position in the CDS
                    row_prediction(row, predictions, 3, minor)
                elif (mutation[-1]=="X"):
                    row=rules.loc[(rules_mutation_type_vector) & (rules.MUTATION=="*X")]
                    #null at any position in the CDS
                    row_prediction(row, predictions, 3, minor)
            else:
                if mutation[-1]=="z":
                    row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="-*z")]
                    #het SNP at any position in the PROM
                    row_prediction(row, predictions, 3, minor)
                elif mutation[-1]=='o':
                    row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="-*o")]
                    #filter fail at any position in the PROM
                    row_prediction(row, predictions, 3, minor)
                elif mutation[-1]=='x':
                    row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="-*x")]
                    #null at any position in the PROM
                    row_prediction(row, predictions, 3, minor)

        # PRIORITY=4: specific mutation at any position in the CDS: only Stop codons make sense for now (e.g. rpoB_*!)
        if mutation[-1] in ['!']:
            if mutation_affects=="CDS":
                row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*!")]
                #stop codon at any position in the CDS
                row_prediction(row, predictions, 4, minor)

        # PRIORITY=5: no change at specific position
        # PRIORITY=6: synonymous mutations at specific location (usually picked up by e.g. L202L, rather than L202=)

        # PRIORITY=7: any het mutation at this specific position in the CDS or PROM  (e.g. rpoB_S450Z or rpoB_c-15z)
        if mutation[-1] in ["Z","z"]:
            row=rules.loc[(rules_mutation_type_vector) & (rules_position_vector) & (rules.MUTATION.str[-1].isin(['Z','z']))]
            #het SNP at specified position in the CDS
            row_prediction(row, predictions, 7, minor)

        # PRIORITY=8: any nonsynoymous mutation at this specific position in the CDS or PROM  (e.g. rpoB_S450? or rpoB_c-15?)
        row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION.str[-1]=="?")]
        #nonsyn SNP at specified position in the CDS
        row_prediction(row, predictions, 8, minor)

    # PRIORITY=9: an exact match
    row=rules.loc[(rules_mutation_type_vector) & (rules.MUTATION==mutation)]
    #exact SNP match
    row_prediction(row, predictions, 9, minor)


def process_indel_variants(mutation_affects: str,
                           predictions: {int: str},
                           rules: pandas.DataFrame,
                           rules_mutation_type_vector: pandas.Series,
                           rules_position_vector: pandas.Series,
                           indel_length: int,
                           indel_type: str,
                           indel_bases: str,
                           position: int,
                           minor: float) -> None:
    '''Apply the cascading rules for INDELs in GARC

    Args:
        mutation_affects (str): Which area this mutation affects. i.e PROM/CDS
        predictions ({int: str}): Dictionary mapping priority -> prediction. This is updated for implict return
        rules (pandas.DataFrame): Rules DataFrame
        rules_mutation_type_vector (pandas.Series): Series for the rules which refer to this mutation type
        rules_position_vector (pandas.Series): Series for the rules which refer to this position
        indel_length (int): Length of the indel
        indel_type (str): Type of the indel. i.e ins/del/fs/indel
        indel_bases (str): Bases inserted/deleted (if applicable)
        position (int): Position of the indel
        minor (float): Float for supporting evidence of minor populations (or None if not a minor population)
    '''
    if mutation_affects == "GENE":
        #Large deletions are priority 1
        large_del(predictions, rules.loc[rules_mutation_type_vector], indel_length, minor)

    # PRIORITY 1: any insertion or deletion in the CDS or PROM (e.g. rpoB_*_indel or rpoB_-*_indel)
    if mutation_affects == "CDS":
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*_indel")]
    else:
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="-*_indel")]
    #any insertion or deletion in the CDS or PROM
    row_prediction(row, predictions, 1, minor)

    # PRIORITY 2: rpoB_*_ins, rpoB_*_del any insertion (or deletion) in the CDS or PROM
    if mutation_affects == "CDS":
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION.isin(["*_ins","*_del"]))]
    else:
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION.isin(["-*_ins","-*_del"]))]
    #any insertion (or deletion) in the CDS or PROM
    row_prediction(row, predictions, 2, minor)

    # PRIORITY 3: any insertion of a specified length (or deletion) in the CDS or PROM (e.g. rpoB_*_ins_2, rpoB_*_del_3, rpoB_-*_ins_1, rpoB_-*_del_200)
    if indel_length is not None and indel_type!="indel":
        if mutation_affects == "CDS":
            row=rules.loc[rules_mutation_type_vector & (rules.MUTATION.isin(["*_"+indel_type+"_"+str(indel_length)]))]
        else:
            row=rules.loc[rules_mutation_type_vector & (rules.MUTATION.isin(["-*_"+indel_type+"_"+str(indel_length)]))]
        #any insertion of a specified length (or deletion) in the CDS or PROM
        row_prediction(row, predictions, 3, minor)

    # PRIORITY=4: any frameshifting insertion or deletion in the CDS (e.g. rpoB_*_fs)
    if indel_length is not None and (indel_length%3)!=0:
        row=rules.loc[rules_mutation_type_vector & (rules.MUTATION=="*_fs")]
        #any frameshifting insertion or deletion in the CDS
        row_prediction(row, predictions, 4, minor)

    # PRIORITY=5: any indel at a specific position in the CDS or PROM (e.g. rpoB_1300_indel or rpoB_-15_indel)
    row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION.str.contains("indel"))]
    #any indel at a specific position in the CDS or PROM (e.g. rpoB_1300_indel or rpoB_-15_indel)
    row_prediction(row, predictions, 5, minor)

    # PRIORITY=6: an insertion (or deletion) at a specific position in the CDS or PROM (e.g. rpoB_1300_ins or rpoB_-15_del)
    if indel_type!="indel":
        row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION.str.contains(indel_type))]
        #any insertion (or deletion) at a specific position in the CDS or PROM (e.g. rpoB_1300_ins or rpoB_-15_del)
        row_prediction(row, predictions, 6, minor)

    # PRIORITY=7: an insertion (or deletion) of a specified length at a specific position in the CDS or PROM (e.g. rpoB_1300_ins_2 or rpoB_-15_del_200)
    if indel_type!="indel" and indel_length is not None:
        row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION==str(position)+"_"+indel_type+"_"+str(indel_length))]
        #an insertion (or deletion) of a specified length at a specific position in the CDS or PROM (e.g. rpoB_1300_ins_2 or rpoB_-15_del_200)
        row_prediction(row, predictions, 7, minor)

    # PRIORITY=8: a frameshifting mutation at a specific position in the CDS (e.g. rpoB_100_fs)
    if indel_length is not None and (indel_length%3)!=0 and position is not None:
        row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION==str(position)+"_fs")]
        #a frameshifting mutation at a specific position in the CDS (e.g. rpoB_100_fs)
        row_prediction(row, predictions, 8, minor)

    # PRIORITY=9: an insertion of a specified series of nucleotides at a position in the CDS or PROM (e.g. rpoB_1300_ins_ca)
    if indel_type=="ins" and indel_length is not None and indel_bases is not None:
        row=rules.loc[rules_mutation_type_vector & rules_position_vector & (rules.MUTATION==str(position)+"_ins_"+indel_bases)]
        #an insertion of a specified series of nucleotides at a position in the CDS or PROM (e.g. rpoB_1300_ins_ca)
        row_prediction(row, predictions, 9, minor)


def parse_mutation(mutation: str) -> (str, str, str, str, int, str, str, str):

    """
    Take a GENE_MUTATION, and determine what type it is, what it affects etc

    Args:
        mutation (str): in the gene_mutation form as defined by GARC e.g. rpoB_S450L or rpoB_1300_ins_3

    Returns:
        position (str): the position of the mutation. This is the amino acid number in a protein, or the number of the nucleotide in the promoter, or a wildcard like '*' or '-*'
        mutation_affects (str): whether this mutation affects the promoter (PROM) or coding sequence (CDS)
        mutation_type (str): is it a SNP or an INDEL?
        indel_type (str): indel, ins, del or fs
        indel_length (int): if sufficient information is given, the number of bases in the INDEL
        indel_bases (str): if sufficient information is given, the bases in the INDEL
        before (str): the REF amino acid
        after (str): the ALT amino acid

    """

    # set default values
    mutation_type=None
    mutation_affects=None
    position=None
    indel_type=None
    indel_length=None
    indel_bases=None
    before=None
    after=None

    # split the mutation on underscore
    cols=mutation.split("_")

    # infer what type it is depending on how many elements it contains
    if len(cols)==1:

        mutation_type="SNP"

        position=int(mutation[1:-1])

        mutation_affects=infer_mutation_affects(position)

        before=mutation[0]
        after=mutation[-1]

    elif len(cols) in [2,3]:
        mutation_type="INDEL"
        #Checking for large deletions
        percentage = False
        if cols[0] == "del":
            try:
                float(cols[1])
                percentage = True
            except ValueError:
                percentage = False

        if percentage:
            mutation_affects = "GENE"
            indel_type = "del"
            indel_length = float(cols[1])
        else:
            try:
                position = int(cols[0])
            except ValueError:
                assert False, "Invalid mutation: " + mutation
            mutation_affects=infer_mutation_affects(position)

            # the third element is one of indel, ins, del or the special case fs
            indel_type=cols[1]

            assert indel_type in ['indel','ins','del','fs'], "form of indel not recognised: "+indel_type

            # if there is a fourth and final element to an INDEL it is either _4 or _ctgc
            if len(cols)==3:
                assert indel_type in ['ins','del'], "form of indel does not make sense when length or bases specified!: "+indel_type
                try:
                    indel_length=int(cols[2])
                except:
                    indel_length=len(cols[2])
                    indel_bases=cols[2]
                    assert 0 not in [c in ['a','t','c','g','z','x'] for c in indel_bases], "only nucleotides of a,t,c,g,z,x are allowed!"

                assert indel_length>0, "makes no sense to have a negative indel length! "+cols[2]

    assert mutation_type in ['INDEL','SNP'], "form of mutation_type not recognised: "+mutation_type

    # insist we've been given an amino acid or a wildcard only
    if mutation_type=="SNP":
        sanity_check_snp(before,after)

    return position, mutation_affects, mutation_type, indel_type, indel_length, indel_bases, before, after

def infer_mutation_affects(position: int) -> str:
    '''Find out which part of the gene this position affects

    Args:
        position (int): Gene position

    Returns:
        str: Either "PROM" or "CDS"
    '''

    if position<0:
        mutation_affects="PROM"
    else:
        mutation_affects="CDS"
    return mutation_affects

def sanity_check_snp(before: str, after: str) -> None:
    '''Sanity check that a given SNP is valid. i.e check that the bases are valid

    Args:
        before (str): Reference bases/AA
        after (str): Alt bases/AA

    Raises:
        AssertationError: Raised in cases in which the SNP is invalid.
    '''

    assert after in ['a','c','t','g','x','z','o','?',"!",'A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','V','W','X','Y','Z'], after+" is not a recognised amino acid or base!"
    assert before in ['a','c','t','g','?','!','A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'], before+" is not a recognised amino acid or base!"

    if before.islower():
        assert after.islower(), "nucleotides must be lowercase!"
        assert before!=after, "makes no sense for the nucleotide to be the same in a mutation!"
    elif before.isupper():
        assert after.isupper() or after=="!", "amino acids must be UPPERCASE!"
