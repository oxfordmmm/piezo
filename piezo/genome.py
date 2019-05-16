#! /usr/bin/env python

import logging, os, gzip, bz2

import pandas, numpy

from Bio import SeqIO
import vcf

class Genome(object):

    def __init__(self,log_file=None,genbank_file=None):

        '''
        Instantiates a genome object by loading a VCF file and storing the whole genome as a numpy array
        '''

        # parse the supplied genbank using BioPython
        GenBankFile=SeqIO.read(genbank_file,"genbank")

        # extract the whole genome sequence (Seq object)
        GenBankFileSeq=GenBankFile.seq

        self.genbank_reference=True

        # store some of the metadata, if there
        self.id=GenBankFile.id
        if 'organism' in GenBankFile.annotations.keys():
            self.organism=GenBankFile.annotations['organism']
        if 'sequence_version' in GenBankFile.annotations.keys():
            self.sequence_version=GenBankFile.annotations['sequence_version']
        if 'source' in GenBankFile.annotations.keys():
            self.source=GenBankFile.annotations['source']
        if 'taxonomy' in GenBankFile.annotations.keys():
            self.taxonomy=GenBankFile.annotations['taxonomy']

        # convert and store it internally as a numpy array of single chars
        self.bases=numpy.array(list(GenBankFileSeq.tomutable()))

        self.bases=numpy.char.lower(self.bases)

        # store how many bases are in the genome
        self.length=len(self.bases)

        # and an array of positions, counting from 1
        self.positions=numpy.arange(0,self.length,1)

        self.genome_string=''.join(self.bases)


    def __repr__(self):

        line=""
        if hasattr(self,'id'):
            line+=self.id+"\n"
        if hasattr(self,'organism'):
            line+=self.organism+"\n"
        if hasattr(self,'sample_name'):
            line+="Sample: "+self.sample_name+"\n"
        else:
            line+="Reference\n"
        line+=str(self.length)+" bases\n"
        line+=self.genome_string[:10]+"...."+self.genome_string[-10:]+"\n"
        return(line)

    def apply_vcf_file(self,vcf_file=None):

        self.genbank_reference=False

        # remember the full path to the VCF file
        self.vcf_file=vcf_file

        # remember the folder path and the name of the passed VCF file
        (self.vcf_folder,self.vcf_filename)=os.path.split(vcf_file)

        filestem, file_extension = os.path.splitext(self.vcf_filename)

        self.sample_name=filestem

        # open the VCF file
        vcf_reader = vcf.Reader(open(self.vcf_file.rstrip(),'r'))

        for record in vcf_reader:

            # find out the position in the genome
            genome_position=int(record.POS)


            # retrieve the details of the row
            row=record.samples[0]

            # find out the reference bases
            ref_bases=record.REF

            # ..and their length
            length_of_variant=len(record.REF)

            # recreate a string of what we expect the reference to be from the genbank file
            # gbk_bases=''.join(reference_genome[genome_position-1:genome_position-1+length_of_variant])
            # assert ref_bases==gbk_bases

            if row.gt_type==2:

                # find out which is the most likely allele in the list
                gt_after=int(row.gt_alleles[1])

                # and hence the alternate bases
                alt_bases=str(record.ALT[gt_after-1])

                # only alter if is a SNP
                if len(ref_bases)==len(alt_bases):

                    for i in alt_bases:

                        self.bases[genome_position]=i.lower()

                        genome_position+=1

        self.genome_string=''.join(self.bases)

    def __sub__(self, other):

        # first store the array of booleans declaring where the arrays are different
        bools_array=self.bases!=other.bases

        ref_array=self.bases[bools_array]

        alt_array=other.bases[bools_array]

        pos_array=self.positions[bools_array]

        result=[]
        for (i,j,k) in zip(pos_array,ref_array,alt_array):
            result.append((i,j,k))

        return(result)

    def save_fasta(self,file_name=None,compression="gzip",compresslevel=2):

        if compression=="gzip":
            OUTPUT=gzip.open(file_name+".gz",'wb',compresslevel=compresslevel)
        elif compression=="bzip2":
            OUTPUT=bz2.open(file_name+".bz2",'wb',compresslevel=compresslevel)
        else:
            OUTPUT=open(file_name,'w')

        header="> "
        if hasattr(self,'id'):
            header+=self.id+"| "
        if hasattr(self,'organism'):
            header+=self.organism+"|  "
        if hasattr(self,'sample_name'):
            header+=self.sample_name
        header+="\n"

        if compression in ['bzip2','gzip']:
            OUTPUT.write(str.encode(header))
            OUTPUT.write(str.encode(self.genome_string))
        else:
            OUTPUT.write(header)
            OUTPUT.write(self.genome_string)

        OUTPUT.close()
