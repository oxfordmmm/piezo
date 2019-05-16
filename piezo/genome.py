#! /usr/bin/env python

import logging, os, gzip, bz2

import pandas, numpy

from Bio import SeqIO
import vcf

class Genome(object):

    def __init__(self,log_file=None,genbank_file=None,fasta_file=None):

        '''
        Instantiates a genome object by loading a VCF file and storing the whole genome as a numpy array
        '''

        assert ((genbank_file is not None) or (fasta_file is not None)), "one of a GenBank file or a FASTA file must be specified!"

        if genbank_file is not None:

            # parse the supplied genbank using BioPython
            GenBankFile=SeqIO.read(genbank_file,"genbank")

            # extract the whole genome sequence (Seq object)
            GenBankFileSeq=GenBankFile.seq

            self.genbank_reference=True
            self.name="Reference"
            self.additional_metadata=None

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

        elif fasta_file is not None:

            if fasta_file.endswith(".gz"):
                INPUT = gzip.open(fasta_file,'rb')
                header=INPUT.readline().decode()
                nucleotide_sequence=INPUT.read().decode()
            elif fasta_file.endswith(".bz2"):
                INPUT = bz2.open(fasta_file,'rb')
                header=INPUT.readline().decode()
                nucleotide_sequence=INPUT.read().decode()
            else:
                INPUT = open(fasta_file,'r')
                header=INPUT.readline()
                nucleotide_sequence=INPUT.read()

            self.genbank_reference=False

            cols=header[1:].split("|")

            if len(cols)>1:
                self.id=cols[0]
                self.organism=cols[1]
                self.name=cols[2]
                if self.name=="Reference":
                    self.genbank_reference=True
            if len(cols)>3:
                self.additional_metadata=cols[3]



            self.bases=numpy.array(list(nucleotide_sequence))

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
        if not self.genbank_reference:
            line+="Sample: "+self.name+"\n"
        else:
            line+="Reference\n"
        if hasattr(self,'path'):
            line+="Path: "+self.path+"\n"
        line+=str(self.length)+" bases\n"
        line+=self.genome_string[:10]+"...."+self.genome_string[-10:]+"\n"
        return(line)

    def apply_vcf_file(self,filename=None):

        self.genbank_reference=False

        # remember the full path to the VCF file
        self.vcf_file=filename

        # remember the folder path and the name of the passed VCF file
        (self.vcf_folder,self.vcf_filename)=os.path.split(self.vcf_file)

        filestem, file_extension = os.path.splitext(self.vcf_filename)

        self.name=filestem

        self.path=self.vcf_folder

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

    def save_fasta(self,filename=None,compression=None,compresslevel=2,additional_metadata=None):

        if compression=="gzip":
            OUTPUT=gzip.open(filename+".gz",'wb',compresslevel=compresslevel)
        elif compression=="bzip2":
            OUTPUT=bz2.open(filename+".bz2",'wb',compresslevel=compresslevel)
        else:
            OUTPUT=open(filename,'w')

        header=">"
        if hasattr(self,'id'):
            header+=self.id+"|"
        if hasattr(self,'organism'):
            header+=self.organism+"|"
        if hasattr(self,'name'):
            header+=self.name
        if additional_metadata is not None:
            header+="|" + additional_metadata
        header+="\n"

        if compression in ['bzip2','gzip']:
            OUTPUT.write(str.encode(header))
            OUTPUT.write(str.encode(self.genome_string))
        else:
            OUTPUT.write(header)
            OUTPUT.write(self.genome_string)

        OUTPUT.close()
