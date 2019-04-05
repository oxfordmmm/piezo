#! /usr/bin/env python

import pkg_resources, os, logging

import numpy, vcf
from Bio import SeqIO

import piezo

class GeneCollection(object):

    """Gene panel class that contains several CRyPTIC gene objects"""

    def __init__(self,species=None,genbank_file=None,log_file=None,gene_panel=None):

        # store the species name, study and instance
        self.species=species
        self.gene_panel=gene_panel

        logging.basicConfig(filename=log_file,level=logging.INFO,format='%(levelname)s, %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

        if self.species=="M. tuberculosis":
            # 4411532 is too small! This assumes no gene name has more than 7 characters
            self.gene_panel_index=numpy.zeros(4420000,dtype='U7')
        else:
            self.gene_panel_index=numpy.zeros(10000000,dtype='U7')

        # remember the config path
        self.config_path = '/'.join(('..','config'))

        # self._parse_genbank_file(pkg_resources.resource_filename("piezo", self.config_path+"/"+genbank_file))
        self._parse_genbank_file(genbank_file)

    def apply_vcf_file(self,vcf_file):

        # remember the full path to the VCF file
        self.vcf_file=vcf_file

        # remember the folder path and the name of the passed VCF file
        (self.vcf_folder,self.vcf_filename)=os.path.split(vcf_file)

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

                # insist that the REF bases indeed match
                vcf_bases=record.REF
                gbk_bases=self.reference_genome[position-1:position-1+len(vcf_bases)].seq
                if vcf_bases!=gbk_bases:
                    logging.warning("REF base(s) mismatch between VCF ("+vcf_bases+") and GENBANK ("+gbk_bases+") at position "+str(position)+" in VCF file "+vcf_file)
                else:

                    # find out what gene we are in
                    gene_name=self.gene_panel_index[position]

                    # be really, really defensive and insist that gt_alleles works as I think it does
                    assert len(row.gt_alleles)==2, "there are more alleles than expected!"

                    (gt_before,gt_after)=(int(row.gt_alleles[0]),int(row.gt_alleles[1]))

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

                    coverage_before=row['COV'][0]
                    coverage_after=row['COV'][gt_after]
                    model_value=row['GT_CONF']

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

                                if before!=after:
                                    self.gene[gene_name].mutate_base(position=position,original_base=before,new_base=after,coverage=[coverage_before,coverage_after],model_score=model_value)

                                # increment the position in the genome
                                position+=1

                        elif len(ref_bases)!=len(alt_bases):

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

                                mutation_name=str(location)+"_indel"

                                self.gene[gene_name].store_indel(mutation=mutation_name,ref=ref_bases,alt=alt_bases,coverage=[coverage_before,coverage_after],model_score=model_value,genome_position=position)


        return(n_hom,n_het,n_ref,n_null)

    def _parse_genbank_file(self,genbank_file):

        self.gene={}

        # read in the M. tuberculosis reference genome
        self.reference_genome=SeqIO.read(genbank_file,'genbank')

        # iterate through all the features in the genomes
        for record in self.reference_genome.features:

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

                    coding_nucleotides=self.reference_genome[start:end]
                    first_nucleotide_index=start+1

                    if strand==1:
                        reverse=False
                        length_promoter=numpy.sum(self.gene_panel_index[start-default_promoter_length:start]=="")
                        self.gene_panel_index[start-length_promoter:end+1]=gene_name
                        promoter_nucleotides=self.reference_genome[start-length_promoter:start]
                    else:
                        reverse=True
                        length_promoter=numpy.sum(self.gene_panel_index[end:end+default_promoter_length]=="")
                        self.gene_panel_index[start:end+1+length_promoter]=gene_name
                        promoter_nucleotides=self.reference_genome[end:end+length_promoter]

                    # store them as a string
                    coding_nucleotides_string=str(coding_nucleotides.seq)
                    promoter_nucleotides_string=str(promoter_nucleotides.seq)

                    # create a Gene object defined using the above information
                    self.gene[gene_name]=piezo.Gene(gene_name=gene_name,coding_nucleotides=coding_nucleotides_string,promoter_nucleotides=promoter_nucleotides_string,first_nucleotide_index=first_nucleotide_index,first_amino_acid_position=first_amino_acid_position,reverse=reverse,element_type=self.gene_panel[gene_name])


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
