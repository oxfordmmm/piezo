#! /usr/bin/env python

import json, re, os, shutil, argparse, logging, pathlib
from datetime import datetime

import pandas
from tqdm import tqdm

import cryptic.genetics
import cryptic.misc

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--database",default='CRyPTIC2',help="the name of the database (default='CRyPTIC2')")
    parser.add_argument("--study",default='V2',help="the name of the study (default='V2')")
    parser.add_argument("--file_type",default="interim",type=str,help="the type of clockwork VCF supplied; either interim or full")
    parser.add_argument("--version",default="0.5.2",type=str,help="the version of clockwork used to create the VCFs")
    parser.add_argument("--testset",action="store_true",help="if specified, will randomly sub-sample the list of identifiers pulled back from CliRes to create a small dataset for testing")
    options = parser.parse_args()

    # be paranoid about what is given to the command line options.
    assert options.file_type in ['interim','ena','full'], 'file_type can only be one of interim or full!'
    assert options.database in ['CRyPTIC1','CRyPTIC2','ENA'], 'only valid entries for database are CRyPTIC1 or CRyPTIC2'
    assert options.study in ['V1','V2','ENA'], 'only valid entries for study is V1 or V2'

    # check the output tree exists (it probably does)
    pathlib.Path('dat/'+options.database+"/"+options.study).mkdir(parents=True, exist_ok=True)

    # construct the path to read the VCFs from based on the command line flags given
    if options.file_type=='interim':
        vcf_path='vcfs/ftp-private.ebi.ac.uk/Interim_VCF_files/Variant_call_version.'+options.version+'/'
    elif options.file_type=='ena':
        vcf_path='vcfs/ftp-private.ebi.ac.uk/ENA_VCF_files/'
    else:
        vcf_path='vcfs/ftp-private.ebi.ac.uk/Internal_releases/1/'

    # check the log folder exists (it probably does)
    pathlib.Path('logs/').mkdir(parents=True, exist_ok=True)

    # open a log file
    datestamp = datetime.strftime(datetime.now(), '%Y-%m-%d_%H%M')
    logging.basicConfig(filename="logs/"+options.database+"-"+options.study+"-vcf-move-"+datestamp+".csv",level=logging.INFO,format='%(levelname)s, %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

    # read in the master JSON file containing the metadata for all the Clockwork VCFs
    with open(vcf_path+'data_by_filename.json') as f:
        Clockwork_VCF_metadata = json.load(f)

    # initialise some global variables
    vcf_counter=0
    categories={}

    # if we are testing, read in the relevant JSON file that the create-phenotype-create-testset.py created
    if options.testset:

        with open("input/"+options.database+"-"+options.study+"-testset-identifiers.json") as INPUT:
            list_of_identifiers=json.load(INPUT)

    # iterate the first time through the Clockwork JSON file to identify any sequencing re-runs
    for vcf_filename in tqdm(Clockwork_VCF_metadata):

        # check that the file actually exists...
        assert os.path.exists(vcf_path+vcf_filename), "VCF file specified in the JSON file does not exist!"

        # if it does, retrieve the key fields for making the UID and path
        (site,subjid,labid,iso,seq_reps)=(Clockwork_VCF_metadata[vcf_filename]['site'],Clockwork_VCF_metadata[vcf_filename]['subject'],Clockwork_VCF_metadata[vcf_filename]['lab_id'],Clockwork_VCF_metadata[vcf_filename]['iso'],Clockwork_VCF_metadata[vcf_filename]['seq_reps'])

        # construct a list (FIXME: could build this off config/CRyPTIC2-V2-path.txt in future)
        list_of_fields=[options.database,options.study,site,subjid,labid,iso,seq_reps]

        # make the UNIQUEID
        uid=cryptic.misc.create_unique_identifier(site,subjid,labid,iso)

        # use the global function to check they are all valid
        (valid,error_text)=cryptic.misc.validate_uid_fields(list_of_fields,data_source='vcf')

        # if not, record to the log file
        if not valid:
            logging.error(uid+error_text)

        # try correcting any problems
        (site,subjid,labid,iso,seq_reps)=cryptic.misc.correct_known_naming_issues(site,subjid.upper(),labid.upper(),iso,seq_reps)

        # construct the same list again, enforcing UPPERCASE
        list_of_fields=[options.database,options.study,site,subjid,labid,iso,seq_reps]

        # re-make the UNIQUEID, correcting any mistakes
        uid=cryptic.misc.create_unique_identifier(site,subjid,labid,iso)

        # use the global function to check if these are all valid
        (valid,error_text)=cryptic.misc.validate_uid_fields(list_of_fields,data_source='vcf')

        # by default, move the VCF file
        proceed=True

        # if there is a problem with the key fields now, don't record to the log file..
        if not valid:
            #  .. instead don't go any further
            proceed=False

        # also stop if we are testing and this UNIQUEID isn't in our testset
        elif options.testset and site+"-"+subjid not in list_of_identifiers.keys():
            proceed=False

        # now check to see if this is the VCF with the longest seq_reps
        # e.g. seqreps.1_2_3 should be used in preference over seqreps.1_2
        if proceed:

            # make the output path if it doesn't already exist
            output_path="dat/"+options.database+"/"+options.study+"/"
            for directory_structure in [site,subjid,labid,iso]:
                output_path+=directory_structure+"/"
                if not os.path.exists(output_path):
                    os.mkdir(output_path)

            # if this is the first time we've encountered this VCF..
            # ..or if the current seqreps is just a longer string than the other VCF file
            if uid not in categories.keys() or len(seq_reps) > len(categories[uid]['SEQREPS']):

                categories[uid]={}
                categories[uid]['OUTPUT_PATH']=output_path
                categories[uid]['STUDYID']=options.database
                categories[uid]['INSTANCE']=options.study
                categories[uid]['UNIQUEID']=uid
                categories[uid]['SITEID']=site
                categories[uid]['SUBJID']=subjid
                categories[uid]['LABID']=labid
                categories[uid]['ISOLATENO']=int(iso)
                categories[uid]['SEQREPS']=seq_reps
                categories[uid]['EBI_ID']=Clockwork_VCF_metadata[vcf_filename]['ebi_internal_data']['id']
                categories[uid]['EBI_ISOLATE_ID']=Clockwork_VCF_metadata[vcf_filename]['ebi_internal_data']['isolate_id']
                categories[uid]['EBI_SAMPLE_ID']=Clockwork_VCF_metadata[vcf_filename]['ebi_internal_data']['sample_id']
                categories[uid]['EBI_PIPELINE_DIR']=Clockwork_VCF_metadata[vcf_filename]['ebi_internal_data']['pipeline_dir']
                categories[uid]['FILENAME']=Clockwork_VCF_metadata[vcf_filename]['filename']
                categories[uid]['VCFSAMPLE']=Clockwork_VCF_metadata[vcf_filename]['vcf_sample']
                categories[uid]['VERSION']=options.version
                categories[uid]['FILETYPE']=options.file_type

    # iterate through the JSON file
    for (uid,categ) in tqdm(categories.items()):

        vcf_filename=categ['FILENAME']
        output_path=categ['OUTPUT_PATH']
        output_filename=uid+"."+options.file_type+".v"+options.version+".vcf"

        try:

            # copy the old file, whilst renaming
            shutil.copy(vcf_path+vcf_filename,output_path+output_filename)

            # instantiate the VCF Treant
            measurement=cryptic.genetics.VCFMeasurement(output_path,new=False,categories=categ,tags=["vcf",'genetics',options.database,options.study,options.file_type,"Clockwork_v"+options.version])

        except:
             logging.error(output_path+output_filename+" - problem copying file into hierarchical structure")

        vcf_counter+=1

    print("A total of %d files have been moved" % vcf_counter)
