#! /usr/bin/env python

import pkg_resources, os, logging, pathlib
from collections import defaultdict

import numpy, pysam
from Bio import SeqIO
from tqdm import tqdm
import piezo
from snpit.snpit import Genotype
from typing import Tuple, List, Dict, IO

class GeneCollection(object):

    """Gene panel class that contains several gene objects"""

    def __init__(self,species=None,genbank_file=None,log_file=None,gene_panel=None,promoter_length=100):

        # store the species name, dictionary of genes and length of promoter
        self.species=species
        self.gene_panel=gene_panel
        self.promoter_length=promoter_length

        # setup a stream to write a log file to
        logging.basicConfig(filename=log_file,level=logging.INFO,format='%(levelname)s, %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

        # load the supplied GenBank file
        self.reference_genome=SeqIO.read(genbank_file,'genbank')

        # convert to a numpy array at the first opportunity since slicing BioPython is between 10 and 50,000 times slower!
        self.reference_genome_sequence=numpy.array([i.lower() for i in str(self.reference_genome.seq)])

        # find out the length of the genome
        self.length_of_genome=len(self.reference_genome_sequence)

        self.gene_panel_index=numpy.zeros(int(1.05*self.length_of_genome),dtype='U10')

        # remember the config path
        self.config_path = '/'.join(('..','config'))

        self._parse_genbank_file()

    def _is_record_invalid(self, ignore_filter: bool, record: pysam.VariantRecord) -> bool:
        return self.gene_panel_index[record.pos]=="" or (
            not ignore_filter and "PASS" not in record.filter.keys()
        )

    @staticmethod
    def get_variant_for_genotype_in_vcf_record(
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
            variant = record.alleles[genotype.call1].lower(),record.alleles[genotype.call2].lower()
        elif genotype.is_null():
            variant = "x"*len(ref_bases),"x"*len(ref_bases)
            print("NULL",ref_bases,variant)
        else:
            raise UnexpectedGenotypeError(
                """Got a genotype for which a Ref/Alt/Null call could not be
                    determined: {}.\nPlease raise this with the developers.""".format(
                    genotype.call()
                )
            )


        return ref_bases,record.pos,variant


    def apply_vcf_file(self,vcf_file,ignore_filter=False, ignore_status=False):

        # remember the full path to the VCF file
        self.vcf_file=vcf_file

        # remember the folder path and the name of the passed VCF file
        (self.vcf_folder,self.vcf_filename)=os.path.split(vcf_file)

        n=defaultdict(int)

        # open the supplied VCF file
        vcf_reader = pysam.VariantFile(self.vcf_file.rstrip())

        # now iterate through the records found in the VCF file
        for record in tqdm(vcf_reader,disable=True):

            # check to see if the position is in one of the genes we are tracking or the filter is not PASS
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
                    genotype = minos_gt_in_wrong_position_fix(record, sample_idx)
                    if genotype is None:
                        raise err

                # if it is simply called reference, note and move on (speed)
                if genotype.is_reference():
                    n['ref']+=1
                    continue

                # find out the variant
                # if it is a het, a tuple is returned for alt_bases
                ref_bases,position,alt_bases = self.get_variant_for_genotype_in_vcf_record(genotype, record)

                # find out what gene we are in
                gene_name=self.gene_panel_index[position]

                if genotype.is_null():
                    print(position,gene_name,ref_bases,alt_bases)

                # insist that the REF bases indeed match
                gbk_bases=''.join(map(str,self.reference_genome_sequence[position-1:position-1+len(ref_bases)]))
                assert ref_bases==gbk_bases, "REF base(s) mismatch between VCF ("+ref_bases+") and GENBANK ("+gbk_bases+") at position "+str(position)+" in VCF file "+vcf_file

                coverage_before=sample_info['COV'][0]

                # store the coverage as a tuple if it is a hetcall
                if genotype.is_heterozygous():
                    coverage_after=(sample_info['COV'][genotype.call1],sample_info['COV'][genotype.call2])
                else:
                    coverage_after=sample_info['COV'][genotype.call1]

                if 'GT_CONF' in sample_info.keys():
                    gt_conf=sample_info['GT_CONF']
                else:
                    gt_conf=None
                if 'GT_CONF_PERCENTILE' in sample_info.keys():
                    gt_conf_percentile=sample_info["GT_CONF_PERCENTILE"]
                else:
                    gt_conf_percentile=None

                if genotype.is_null():
                    n['null']+=1
                    print("NULL CALL:", ref_bases, alt_bases, coverage_before, coverage_after,gt_conf,gt_conf_percentile)
                elif genotype.is_alt():
                    n['hom']+=1
                    # print("HOM CALL:", ref_bases, alt_bases, coverage_before, coverage_after,gt_conf,gt_conf_percentile)
                elif genotype.is_heterozygous():
                    n['het']+=1
                    if ref_bases!=alt_bases[0]:
                        print("HET CALL:", ref_bases, alt_bases, coverage_before, coverage_after,gt_conf,gt_conf_percentile)
                    # FIXME:

                # print("COVERAGE",coverage_before,coverage_after)

                # automatically not true if a het call since alt_bases will be a tuple
                # true for a SNP
                if len(ref_bases)==len(alt_bases):

                    # most of the time these will be SNPs so the loop will iterate just once
                    for before,after in zip(ref_bases,alt_bases):

                        # the mutate base is setup to handle one base at a time, so call it each time around the loop
                        # note that this will be called even if before==after, but that is needed to deal with k-mers where not all bases have been mutated
                        # the mutate_base() method does assert that the original_base matches what is in the appropriate GenBank file

                        # if before!=after:

                            # print(gene_name,position,before,after,coverage_before,coverage_after,gt_conf,ref_bases,alt_bases)
                            # try:
                            # self.gene[gene_name].mutate_base(position=position,original_base=before,new_base=after,coverage=[coverage_before,coverage_after],model_score=gt_conf)
                            # except:
                            #     print(gene_name,before,after)

                        # increment the position in the genome
                        position+=1

                #
                # if not genotype.is_heterozygous():
                #
                #     # alt_bases
                #     pass
                #
                # if row.gt_type==2 or row['GT']=='./.' or row.gt_type==1:
                #
                #     # find out what the bases are in tge GenBank reference genome
                #     ref_bases=record.REF
                #
                #     # insert 'x' base for low-coverage nulls
                #     if row['GT']=='./.':
                #         alt_bases='x' * len(ref_bases)
                #
                #     # insert 'z' base for hets
                #     elif row.gt_type==1:
                #         alt_bases='z' * len(ref_bases)
                #
                #     # otherwise it must be a homozygous call so use the genotype to determine the correct alt
                #     else:
                #         alt_bases=str(record.ALT[gt_after-1])
                #
                #     # we are now in a position to either make a mutation (SNP) or remember an INDEL
                #     # the majority of the below will be SNPs, but some will be e.g. 5-mers in gyrA with 2 SNPs in close proximity
                #     if len(ref_bases)==len(alt_bases):
                #
                #         # most of the time these will be SNPs so the loop will iterate just once
                #         for before,after in zip(ref_bases,alt_bases):
                #
                #             # the mutate base is setup to handle one base at a time, so call it each time around the loop
                #             # note that this will be called even if before==after, but that is needed to deal with k-mers where not all bases have been mutated
                #             # the mutate_base() method does assert that the original_base matches what is in the appropriate GenBank file
                #
                #             if before!=after:
                #                 try:
                #                     self.gene[gene_name].mutate_base(position=position,original_base=before,new_base=after,coverage=[coverage_before,coverage_after],model_score=model_value)
                #                 except:
                #                     print(gene_name)
                #
                #             # increment the position in the genome
                #             position+=1
                #
                #     elif len(ref_bases)!=len(alt_bases):
                #
                #         # find out where the mutation is
                #         (type,location)=self.gene[gene_name].return_location(position)
                #
                #         # be defensive; the above only returns None if the position isn't in the gene
                #         if type!=None:
                #
                #             element_type=self.gene[gene_name].element_type
                #
                #             if type=="CDS":
                #                 cds=True
                #                 promoter=False
                #             elif type=="PROMOTER":
                #                 cds=False
                #                 promoter=True
                #             else:
                #                 raise ValueError("Only CDS or PROMOTER can be returned: was given "+type)
                #
                #             mutation_name=str(location)+"_indel"
                #
                #             self.gene[gene_name].store_indel(mutation=mutation_name,ref=ref_bases,alt=alt_bases,coverage=[coverage_before,coverage_after],model_score=model_value,genome_position=position)

        # force all the gene strings etc to be rebuilt now that the mutations have been applied
        self._update_genes()
        print(n['hom'],n['het'],n['ref'],n['null'])
        return(n['hom'],n['het'],n['ref'],n['null'])

    def _update_genes(self):

        for gene_name in self.gene_panel:
            self.gene[gene_name].update()

    def _parse_genbank_file(self):

        self.gene={}

        # iterate through all the features in the genomes
        for record in self.reference_genome.features:

            # check that the record is a Coding Sequence and it is also a gene
            if record.type not in ['CDS','rRNA']:
                continue

            # if it is a gene
            if 'gene' in record.qualifiers.keys():
                if any(i in record.qualifiers['gene'] for i in self.gene_panel):
                    gene_name=record.qualifiers['gene'][0]
                else:
                    continue

            elif 'locus_tag' in record.qualifiers.keys():
                if any(i in record.qualifiers['locus_tag'] for i in self.gene_panel):
                    gene_name=record.qualifiers['locus_tag'][0]
                else:
                    continue

            # the start and end positions in the reference genome
            start=record.location.start.position
            end=record.location.end.position

            # and finally the direction of the strand (which if -1 means reverse complement)
            strand=record.location.strand.real

            # retrieve the number of the first protein coding codon (usually 1)
            if self.gene_panel[gene_name] in ['GENE','LOCUS']:
                codon_start=record.qualifiers['codon_start']
                first_amino_acid_position=int(codon_start[0])
                promoter_length=self.promoter_length
            elif self.gene_panel[gene_name]=='RNA':
                first_amino_acid_position=1
                promoter_length=0
            else:
                raise Exception("gene does not have correct type specified in gene panel file")

            coding_nucleotides=self.reference_genome_sequence[start:end]
            first_nucleotide_index=start+1

            if strand==1:
                reverse=False
                length_promoter=numpy.sum(self.gene_panel_index[start-promoter_length:start]=="")
                self.gene_panel_index[start-length_promoter:end+1]=gene_name
                promoter_nucleotides=self.reference_genome_sequence[start-length_promoter:start]
            else:
                reverse=True
                length_promoter=numpy.sum(self.gene_panel_index[end:end+promoter_length]=="")
                self.gene_panel_index[start:end+1+length_promoter]=gene_name
                promoter_nucleotides=self.reference_genome_sequence[end:end+length_promoter]

            # convert them both into strings
            coding_nucleotides_string="".join(map(str,coding_nucleotides))
            promoter_nucleotides_string="".join(map(str,promoter_nucleotides))

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

def minos_gt_in_wrong_position_fix(record, sample_idx):
    """A version of minos had GT in the second column instead of the first"""
    info = str(record).strip().split("\t")[9 + sample_idx]
    for field in info.split(":"):
        if "/" in field:
            return Genotype.from_string(field)
