#! /usr/bin/env python

import pkg_resources, os, logging, pathlib

import numpy, pysam
from Bio import SeqIO
from tqdm import tqdm
import piezo
from snpit.genotype import Genotype
from typing import Tuple, List, Dict, IO

class GeneCollection(object):

    """Gene panel class that contains several gene objects"""

    def __init__(self,species=None,genbank_file=None,log_file=None,gene_panel=None,promoter_length=100):

        # store the species name, dictionary of genes and flags to control VCF reading behaviour
        self.species=species
        self.gene_panel=gene_panel
        self.promoter_length=promoter_length

        logging.basicConfig(filename=log_file,level=logging.INFO,format='%(levelname)s, %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

        # load the supplied GenBank file
        self.reference_genome=SeqIO.read(genbank_file,'genbank')

        # find out the length of the genome
        self.length_of_genome=len(self.reference_genome.seq)

        # 4411532 is too small! This assumes no gene name has more than 9 characters
        self.gene_panel_index=numpy.zeros(int(self.length_of_genome*1.05),dtype='U9')

        self._parse_genbank_file(genbank_file)

    def _is_record_invalid(self, gene, ignore_filter, record: pysam.VariantRecord) -> bool:
        return gene=="" or (
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
        if genotype.is_reference():
            variant = ""
        elif genotype.is_heterozygous():
            variant = 'z'
        elif genotype.is_alt():
            alt_call = max(genotype.call())
            variant = record.alleles[alt_call]
        elif genotype.is_null():
            variant = "-"
        else:
            raise UnexpectedGenotypeError(
                """Got a genotype for which a Ref/Alt/Null call could not be
                    determined: {}.\nPlease raise this with the developers.""".format(
                    genotype.call()
                )
            )
        return variant


    def apply_vcf_file(self,vcf_file,ignore_filter=False, ignore_status=False):

        # remember the full path to the VCF file
        self.vcf_file=vcf_file

        # remember the folder path and the name of the passed VCF file
        (self.vcf_folder,self.vcf_filename)=os.path.split(vcf_file)

        n_null=0
        n_het=0
        n_ref=0
        n_hom=0

        # open the VCF file from the EBI
        vcf_reader = pysam.VariantFile(self.vcf_file.rstrip())

        MUTATIONS_dict={}
        MUTATIONS_counter=0
        EFFECTS_dict={}
        EFFECTS_counter=0

        # now iterate through the records found in the VCF file
        for record in vcf_reader:

            putative_gene=self.gene_panel_index[record.pos]

            # check to see if the position is in one of the genes we are tracking or the filter is not PASS
            if self._is_record_invalid(putative_gene,ignore_filter,record):
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

                variant = self.get_variant_for_genotype_in_vcf_record(genotype, record)

                if not variant:
                    continue

                if genotype.is_null():
                    n_null+=1
                elif genotype.is_reference():
                    n_ref+=1
                elif genotype.is_alt():
                    n_hom+=1
                elif genotype.is_heterozygous():
                    n_het+=1

                # find out our reference genome position
                position=int(record.pos)

                # insist that the REF bases indeed match
                vcf_bases=record.ref
                gbk_bases=self.reference_genome[position-1:position-1+len(vcf_bases)].seq

                if vcf_bases!=gbk_bases:
                    logging.warning("REF base(s) mismatch between VCF ("+vcf_bases+") and GENBANK ("+gbk_bases+") at position "+str(position)+" in VCF file "+vcf_file)
                    continue
                print(genotype.call1, genotype.call2)
                print(sample_info.keys())
                # find out what gene we are in
                # gene_name=

                print(sample_info['COV'])
                print(sample_info['DP'])
                print(sample_info['GT_CONF_PERCENTILE'])

                coverage_before=sample_info['COV'][0]
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
                                try:
                                    self.gene[gene_name].mutate_base(position=position,original_base=before,new_base=after,coverage=[coverage_before,coverage_after],model_score=model_value)
                                except:
                                    print(gene_name)


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

        # iterate through all the features in the genomes
        for record in self.reference_genome.features:

            # check that the record is a Coding Sequence and it is also a gene
            if record.type in ['CDS','rRNA']:

                found_record=False

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

                coding_nucleotides=self.reference_genome[start:end]
                first_nucleotide_index=start+1

                if strand==1:
                    reverse=False
                    length_promoter=numpy.sum(self.gene_panel_index[start-promoter_length:start]=="")
                    self.gene_panel_index[start-length_promoter:end+1]=gene_name
                    promoter_nucleotides=self.reference_genome[start-length_promoter:start]
                else:
                    reverse=True
                    length_promoter=numpy.sum(self.gene_panel_index[end:end+promoter_length]=="")
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

def minos_gt_in_wrong_position_fix(record, sample_idx):
    """A version of minos had GT in the second column instead of the first"""
    info = str(record).strip().split("\t")[9 + sample_idx]
    for field in info.split(":"):
        if "/" in field:
            return Genotype.from_string(field)
