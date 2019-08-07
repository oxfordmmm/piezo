import gzip, os, pickle

import numpy

import pysam
from Bio import SeqIO
from tqdm import tqdm

from piezo import Genotype

class Genome(object):

    def __init__(self,genbank_file=None,fasta_file=None):

        '''
        Instantiates a genome object by loading a VCF file and storing the whole genome as a numpy array

        Args:
            genbank_file (str): path to the GenBank file to build the reference genome
            fasta_file (str): path to the FASTA file to build the reference genome
        '''

        assert ((genbank_file is not None) or (fasta_file is not None)), "one of a GenBank file or a FASTA file must be specified!"

        self.id=""
        self.organism=""
        self.sample_name=""
        self.additional_metadata=None

        # load the specified GenBank file
        if genbank_file is not None:

            # create the genbank file and store in a BioPython object
            reference_genome=SeqIO.read(genbank_file,'genbank')

            self.sample_name="Reference Genome"

            # convert to a numpy array at the first opportunity since slicing BioPython is between 10 and 50,000 times slower!
            self.sequence=numpy.array([i.lower() for i in str(reference_genome.seq)])

            # store some of the metadata, if it is present
            self.id=reference_genome.id

            if 'organism' in reference_genome.annotations.keys():
                self.organism=reference_genome.annotations['organism']
            if 'sequence_version' in reference_genome.annotations.keys():
                self.sequence_version=reference_genome.annotations['sequence_version']
            if 'source' in reference_genome.annotations.keys():
                self.source=reference_genome.annotations['source']
            if 'taxonomy' in reference_genome.annotations.keys():
                self.taxonomy=reference_genome.annotations['taxonomy']

        # otherwise there must be a FASTA file so load that instead
        elif fasta_file is not None:

            header,nucleotide_sequence=self._load_fastafile(fasta_file)

            cols=header[1:].split("|")
            if len(cols)>1:
                self.id=cols[0]
                self.organism=cols[1]
                self.sample_name=cols[2]
            if len(cols)>3:
                self.additional_metadata=cols[3]

            self.sequence=numpy.array(list(nucleotide_sequence))

        # insist that bases are lower case
        self.sequence=numpy.char.lower(self.sequence)

        self.bases_integer_lookup, self.integers = numpy.unique(self.sequence, return_inverse=True)

        # store the length of the genome
        self.length=len(self.sequence)

        # create a diploid version of the sequence to handle HET calls later
        self.diploid=(self.sequence,self.sequence)

        # create an array of the genome indices
        self.index=numpy.arange(1,self.length+1)

    def __repr__(self):

        '''
        Overload the print function to write a summary of the genome.
        '''

        output=""
        if hasattr(self,'id'):
            output+=self.id+"\n"
        if hasattr(self,'organism'):
            output+=self.organism+"\n"
        if hasattr(self,'sample_name'):
            output+=self.sample_name+"\n"
        output+=str(self.length)+" bases\n"
        output+=''.join(i for i in self.sequence[0:3])
        output+="..."
        output+=''.join(i for i in self.sequence[-3:])

        return(output)

    def __sub__(self,other):

        """
        Overload the subtraction operator so it returns a tuple of the differences between the two genomes
        """

        assert self.length==other.length, "genomes must have the same length!"

        mask=self.sequence!=other.sequence

        return(self.sequence[mask],self.index[mask],other.sequence[mask])

    def calculate_snp_distance(self,other):
        return (numpy.count_nonzero(self.sequence!=other.sequence))

    def apply_vcf_file(self,vcf_file=None,ignore_filter=False, ignore_status=False,show_progress_bar=False):
        """
        Load a VCF file and apply the variants to the whole genome sequence.

        Args:
            vcf_file (str): path to the VCF file to be loaded and applied
            ignore_filter (bool): whether to ignore the FILTER column in the VCF file (Clockwork hasn't always written it correctly)
            ignore_status (bool): ditto
            show_progress_bar (bool): whether to draw a nice tqdm progress bar (False by default)
        """

        # since we are now applying a VCF file, it makes sense to create these numpy arrays
        self.coverage=numpy.zeros(self.length)
        self.model_score=numpy.zeros(self.length)
        self.model_percentile=numpy.zeros(self.length)

        # split and remember the path, filename and stem of the VCF file
        (self.vcf_folder,self.vcf_file_name)=os.path.split(vcf_file)
        self.vcf_file_stem, file_extension = os.path.splitext(self.vcf_file_name)

        # assume that the sample name is the filestem and remember
        self.sample_name=self.vcf_file_stem

        # open the supplied VCF file
        # note that this will read in bgzip compressed vcf files (from htslib) but not gzip compressed files
        # even though the file extension is the same
        vcf_reader = pysam.VariantFile(vcf_file.rstrip())

        # now iterate through the records found in the VCF file
        for record in tqdm(vcf_reader,disable=not(show_progress_bar)):

            # check to see the filter is not PASS
            if self._is_record_invalid(ignore_filter,record):
                continue

            for sample_idx, (sample_name, sample_info) in enumerate(
                record.samples.items()
            ):

                if not ignore_status and sample_info["STATUS"] == "FAIL":
                    continue

                try:
                    genotype = Genotype(*sample_info["GT"])
                except TypeError as err:
                    genotype = self._minos_gt_in_wrong_position_fix(record, sample_idx)
                    if genotype is None:
                        raise err

                ref_bases,position,alt_bases = self._get_variant_for_genotype_in_vcf_record(genotype, record)

                if alt_bases=="":
                    continue

                if len(ref_bases)==len(alt_bases):

                    coverage=sample_info['COV'][genotype.call1]
                    model_score=sample_info['GT_CONF']
                    if 'GT_CONF_PERCENTILE' in sample_info.keys():
                        model_percentile=sample_info['GT_CONF_PERCENTILE']

                    for before,after in zip(ref_bases,alt_bases):

                        if before!=after:
                            self.sequence[position-1]=after
                            self.coverage[position-1]=coverage
                            self.model_score[position-1]=model_score
                            if 'GT_CONF_PERCENTILE' in sample_info.keys():
                                self.model_percentile[position-1]=model_percentile

                        # increment the position in the genome
                        position+=1

    def save_fasta(self,filename=None,compression=False,compresslevel=2,additional_metadata=None,chars_per_line=70,nucleotides_uppercase=True):

        '''
        Save the genome as a FASTA file.

        Args:
            filename (str): path of the output file
            compression (bool): If True, save compressed using gzip. (bzip2 is too slow)
            compresslevel (0-9): the higher the number, the harder the algorithm tries to compress but it takes longer. Default is 2.
            additional_metadata (str): will be added to the header of the FASTA file
            chars_per_line (int): the number of characters per line. Default=70. Must be either a positive integer or None (i.e. no CRs)
        '''

        # check the arguments are well formed
        assert compression in [True,False]
        assert nucleotides_uppercase in [True,False]
        assert compresslevel in range(1,10), "compresslevel must be in range 1-9!"
        if chars_per_line is not None:
            assert chars_per_line > 0, "number of characters per line in the FASTA file must be an integer!"

        # check the specified fileextension to see if the FASTA file needs compressing
        if compression:
            OUTPUT=gzip.open(filename+".gz",'wb',compresslevel=compresslevel)
        else:
            OUTPUT=open(filename,'w')

        # create the header line for the FASTA file using "|" as delimiters
        header=">"
        if hasattr(self,'id'):
            header+=self.id+"|"
        if hasattr(self,'organism'):
            header+=self.organism+"|"
        header+=self.sample_name
        if additional_metadata is not None:
            header+="|" + additional_metadata
        header+="\n"

        # create a string of the genome
        genome_string=''.join(self.sequence)

        # insert carriage returns so it looks pretty in the file...
        output_string=self._insert_newlines(genome_string,every=chars_per_line)
        output_string+="\n"

        # set the case accordingly
        if nucleotides_uppercase:
            output_string=output_string.upper()
        else:
            output_string=output_string.lower()

        # write out the FASTA files
        if compression:
            OUTPUT.write(str.encode(header))
            OUTPUT.write(str.encode(output_string))
        else:
            OUTPUT.write(header)
            OUTPUT.write(output_string)

        OUTPUT.close()

    @staticmethod
    def _get_variant_for_genotype_in_vcf_record(
        genotype: Genotype, record: pysam.VariantRecord
    ) -> str:
        """Retrieves the variant a genotype maps to for a given record.

        Args:
            genotype: The genotype call for the sample.
            record: A VCF record object.
        Returns:
            str: A hyphen if the call is null (ie ./.) or the alt variant if
            the call is alt. Returns an empty string if the call is ref or heterozygous.
        """

        #  find out what the reference bases are
        ref_bases=record.ref.lower()

        if genotype.is_reference():
            variant = ""
        elif genotype.is_heterozygous():
            variant = record.alleles[genotype.call1].lower(),record.alleles[genotype.call2].lower()
        elif genotype.is_alt():
            variant = record.alleles[genotype.call1].lower()
        elif genotype.is_null():
            variant = "x"*len(ref_bases)
        else:
            raise UnexpectedGenotypeError(
                """Got a genotype for which a Ref/Alt/Null call could not be
                    determined: {}.\nPlease raise this with the developers.""".format(
                    genotype.call()
                )
            )


        return ref_bases,record.pos,variant

    def _is_record_invalid(self, ignore_filter: bool, record: pysam.VariantRecord) -> bool:
        """
        Simple private method for parsing VCF record
        """
        return ( not ignore_filter and "PASS" not in record.filter.keys() )

    def _minos_gt_in_wrong_position_fix(self,record, sample_idx):

        """
        Hacky private method to fix a minos mistake

        (A version of minos had GT in the second column instead of the first)
        """

        info = str(record).strip().split("\t")[9 + sample_idx]
        for field in info.split(":"):
            if "/" in field:
                return Genotype.from_string(field)


    @staticmethod
    def _load_fastafile(fasta_file):
        """
        Loads the fasta file whether uncompressed or compressed by gzip

        Args:
            fasta_file(path): the path to the fasta file

        Returns:
            header (str): the first line of the FASTA file which will hopefully contain some metadata
            nucleotide_sequence (str): the nucleotide sequence as a string
        """

        # check if it is compressed and load it accordingly
        if fasta_file.endswith(".gz"):
            INPUT = gzip.open(fasta_file,'rb')
            header=INPUT.readline().decode()
            nucleotide_sequence=INPUT.read().decode()
        else:
            INPUT = open(fasta_file,'r')
            header=INPUT.readline()
            nucleotide_sequence=INPUT.read()

        nucleotide_sequence=nucleotide_sequence.replace('\n','')

        return(header,nucleotide_sequence)

    @staticmethod
    def _insert_newlines(string, every=70):
        '''
        Simple private method for inserting a carriage return every N characters into a long string.

        Args:
            string (str): the string to insert carriage returns
            every (int): how many characters between each carriage return
        '''

        assert every>0, "every must be an integer greater than zero"

        assert len(string)>1, "string is too short!"

        return '\n'.join(string[i:i+every] for i in range(0, len(string), every))
