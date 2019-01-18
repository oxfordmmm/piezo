#! /usr/bin/env python

import pkg_resources, os, logging, pathlib
from datetime import datetime

import numpy, vcf
from Bio import SeqIO

from snpit import snpit

import cryptic.genetics

class GeneCollection(object):

    """Gene panel that contains several CRyPTIC gene objects"""

    def __init__(self,species=None,genbank_file=None,database=None,study=None,gene_panel=None):

        # store the species name, study and instance
        self.species=species
        self.database=database
        self.study=study
        self.instance=study
        self.gene_panel=gene_panel

        # check the log folder exists (it probably does)
        pathlib.Path('logs/').mkdir(parents=True, exist_ok=True)

        datestamp = datetime.strftime(datetime.now(), '%Y-%m-%d_%H%M')
        logging.basicConfig(filename="logs/"+self.database+"-"+self.study+"-cryptic-gene-collection-"+datestamp+".csv",level=logging.INFO,format='%(levelname)s, %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

        if self.species=="M. tuberculosis":
            # 4411532 is too small! This assumes no gene name has more than 7 characters
            self.gene_panel_index=numpy.zeros(4420000,dtype='U7')
        else:
            self.gene_panel_index=numpy.zeros(10000000,dtype='U7')

        # remember the config path
        self.config_path = '/'.join(('..','config'))

        # load and store the gene panel file in a dictionary
        # self.gene_panel=self._read_genepanel_file(pkg_resources.resource_filename("cryptic", self.config_path+"/"+self.database+"-"+self.study+"-gene_panel.csv"))

        # print((self.gene_panel))
        # for i in self.gene_panel:
        #     print(i, self.gene_panel[i])

        self._parse_genbank_file(pkg_resources.resource_filename("cryptic", self.config_path+"/"+genbank_file))

    def apply_vcf_file(self,vcf_file):

        # remember the full path to the VCF file
        self.vcf_file=vcf_file

        # remember the folder path and the name of the passed VCF file
        (self.vcf_folder,self.vcf_filename)=os.path.split(vcf_file)

        (self.UniqueID,self.site_id,self.subject_id,self.lab_id,self.iso_id) = self._parse_vcf_filename(self.vcf_filename)

        # find out what the species, lineage etc are
        (self.species,self.lineage,self.sublineage,self.lineage_percentage)=self._determine_lineage()

        n_null=0
        n_het=0
        n_ref=0
        n_hom=0

        # open the VCF file from the EBI
        vcf_reader = vcf.Reader(open(self.vcf_file.rstrip(),'r'))

        MUTATIONS_dict={}
        MUTATIONS_counter=0
        EFFECTS_dict={}
        EFFECTS_counter=0

        # now iterate through the records found in the VCF file
        for record in vcf_reader:

            # be defensive and check we've only got one row per sample as expected
            assert len(record.samples)==1, "Row in VCF with more than one sample"

            # if so, pull out the row
            row=record.samples[0]

            # record some basic statistics
            if row['GT']=='./.':
                n_null+=1
            elif row.gt_type==0:
                n_ref+=1
            elif row.gt_type==1:
                n_het+=1
            elif row.gt_type==2:
                n_hom+=1

            # find out our reference genome position
            position=int(record.POS)

            # and continue only if the mutation occurs in the one of the genes from the panel and if the call was actually made
            if self.gene_panel_index[position]!="" and row.called:

                # find out what gene we are in
                gene_name=self.gene_panel_index[position]

                # be really, really defensive and insist that gt_nums and gt_type works as I think it does
                gt_before=int(row.gt_nums.split("/")[0])
                gt_after=int(row.gt_nums.split("/")[1])

                # also be defensive about the relatioship between gt_type and gt_before/after
                # these should all be 0/0 as are ref calls i.e. leave as reference
                if row.gt_type==0:
                    assert gt_before==gt_after==0, "ref calls not working as expected in VCF file"

                # these are het calls
                elif row.gt_type==1:
                    assert gt_before!=gt_after, "het calls not working as expected in VCF file"

                # whilst these are the homozygous variants
                elif row.gt_type==2:
                    assert gt_before==gt_after

                else:
                    raise ValueError("gt_type is not one of 0,1 or 2 but is instead "+str(row.gt_type))

                # ignore ref calls (0/0) and only consider homozygous, null and hets
                if row.gt_type==2 or row['GT']=='./.' or row.gt_type==1:

                    # find out what the bases are in tge GenBank reference genome
                    ref_bases=record.REF

                    # insert 'x' base for low-coverage nulls
                    if row['GT']=='./.':
                        alt_bases='x' * len(ref_bases)

                    # insert 'z' base for hets
                    elif row.gt_type==1:
                        alt_bases='z' * len(ref_bases)

                    # otherwise it must be a homozygous call so use the genotype to determine the correct alt
                    else:
                        alt_bases=str(record.ALT[gt_after-1])

                    # we are now in a position to either make a mutation (SNP) or remember an INDEL
                    # the majority of the below will be SNPs, but some will be e.g. 5-mers in gyrA with 2 SNPs in close proximity
                    if len(ref_bases)==len(alt_bases):

                        # most of the time these will be SNPs so the loop will iterate just once
                        for before,after in zip(ref_bases,alt_bases):

                            # the mutate base is setup to handle one base at a time, so call it each time around the loop
                            # note that this will be called even if before==after, but that is needed to deal with k-mers where not all bases have been mutated
                            # the mutate_base() method does assert that the original_base matches what is in the appropriate GenBank file
                            self.gene[gene_name].mutate_base(position=position,original_base=before,new_base=after)

                            # increment the position in the genome
                            position+=1

                    elif len(ref_bases)!=len(alt_bases):

                        # if an insertion, this will be positive. A deletion will be negative
                        length_of_indel=len(alt_bases)-len(ref_bases)

                        # find out where the mutation is
                        (type,location)=self.gene[gene_name].return_location(position)

                        # be defensive; the above only returns None if the position isn't in the gene
                        if type!=None:

                            element_type=self.gene[gene_name].element_type

                            if type=="CDS":
                                cds=True
                                promoter=False
                            elif type=="PROMOTER":
                                cds=False
                                promoter=True
                            else:
                                raise ValueError("Only CDS or PROMOTER can be returned: was given "+type)

                            if length_of_indel>0:
                                insertion=True
                                deletion=False
                                mutation_name=str(location)+"_ins_"+str(abs(length_of_indel))
                            else:
                                insertion=False
                                deletion=True
                                mutation_name=str(location)+"_del_"+str(abs(length_of_indel))

                            self.gene[gene_name].store_indel(mutation=mutation_name,ref=ref_bases,alt=alt_bases)


        return(n_hom,n_het,n_ref,n_null)


    # def _read_genepanel_file(self,gene_panel_file):
    #
    #     gene_panel={}
    #     with open(gene_panel_file) as file:
    #         # read and discard the first line as the header
    #         file.readline()
    #         for line in file:
    #             cols=line.split(',')
    #             gene_panel[cols[0]]=cols[1].rstrip().upper()
    #
    #     return(gene_panel)

    def _parse_genbank_file(self,genbank_file):

        self.gene={}

        # read in the M. tuberculosis reference genome
        reference_genome=SeqIO.read(genbank_file,'genbank')

        # iterate through all the features in the genomes
        for record in reference_genome.features:

            # check that the record is a Coding Sequence and it is also a gene
            if record.type in ['CDS','rRNA']:

                found_record=False

                # if it is a gene
                if 'gene' in record.qualifiers.keys():
                    if any(i in record.qualifiers['gene'] for i in self.gene_panel):
                        gene_name=record.qualifiers['gene'][0]
                        found_record=True

                elif 'locus_tag' in record.qualifiers.keys():
                    if any(i in record.qualifiers['locus_tag'] for i in self.gene_panel):
                        gene_name=record.qualifiers['locus_tag'][0]
                        found_record=True

                if found_record:

                    # the start and end positions in the reference genome
                    start=record.location.start.position
                    end=record.location.end.position

                    # and finally whether the strand (which if -1 means reverse complement)
                    strand=record.location.strand.real

                    # retrieve the number of the first protein coding codon (usually 1)
                    if self.gene_panel[gene_name] in ['GENE','LOCUS']:
                        codon_start=record.qualifiers['codon_start']
                        first_amino_acid_position=int(codon_start[0])
                        default_promoter_length=100
                    elif self.gene_panel[gene_name]=='RNA':
                        first_amino_acid_position=1
                        default_promoter_length=0
                    else:
                        raise Exception("gene does not have correct type specified in gene panel file")

                    coding_nucleotides=reference_genome[start:end]
                    first_nucleotide_index=start+1

                    if strand==1:
                        reverse=False
                        length_promoter=numpy.sum(self.gene_panel_index[start-default_promoter_length:start]=="")
                        self.gene_panel_index[start-length_promoter:end+1]=gene_name
                        promoter_nucleotides=reference_genome[start-length_promoter:start]
                    else:
                        reverse=True
                        length_promoter=numpy.sum(self.gene_panel_index[end:end+default_promoter_length]=="")
                        self.gene_panel_index[start:end+1+length_promoter]=gene_name
                        promoter_nucleotides=reference_genome[end:end+length_promoter]

                    # store them as a string
                    coding_nucleotides_string=str(coding_nucleotides.seq)
                    promoter_nucleotides_string=str(promoter_nucleotides.seq)

                    # create a Gene object defined using the above information
                    self.gene[gene_name]=cryptic.genetics.Gene(gene_name=gene_name,coding_nucleotides=coding_nucleotides_string,promoter_nucleotides=promoter_nucleotides_string,first_nucleotide_index=first_nucleotide_index,first_amino_acid_position=first_amino_acid_position,reverse=reverse,element_type=self.gene_panel[gene_name])


    def _determine_lineage(self):

        tb=snpit(input_file=self.vcf_file)

        return (tb.species,tb.lineage,tb.sublineage,tb.percentage)

    def _parse_vcf_filename(self,filename):

        # parse the filename
        cols=filename.split('.')
        site_id=cols[1]
        subject_id=cols[3]
        lab_id=cols[5]
        iso_id=int(cols[7])
        seq_reps=cols[9]

        UniqueID="site."+site_id+".subj."+subject_id+".lab."+lab_id+".iso."+str(iso_id) #+".seq_reps."+seq_reps

        return(UniqueID,site_id,subject_id,lab_id,iso_id)
