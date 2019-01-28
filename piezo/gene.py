#! /usr/bin/env python

from datetime import datetime
import logging

import numpy

# FIXME: problems with rrs, mfpB

class Gene(object):

    """Gene object that uses underlying numpy arrays"""

    def __init__(self,gene_name=None,coding_nucleotides=None,promoter_nucleotides=None,first_amino_acid_position=1,first_nucleotide_index=None,reverse=False,element_type=None):

        self.gene_name=gene_name
        self.element_type=element_type
        self.coding_nucleotides_string=coding_nucleotides.lower()
        self.promoter_nucleotides_string=promoter_nucleotides.lower()
        self.first_nucleotide_index=first_nucleotide_index
        self.reverse=reverse

        # create an empty dict to store indels in later
        self.indels={}

        datestamp = datetime.strftime(datetime.now(), '%Y-%m-%d_%H%M')
        logging.basicConfig(filename="logs/piezo-gene-"+gene_name+"-"+datestamp+".csv",level=logging.INFO,format='%(levelname)s, %(message)s', datefmt='%a %d %b %Y %H:%M:%S')

        # how many nucleotides have we been given?
        self.number_coding_nucleotides=len(self.coding_nucleotides_string)
        self.last_nucleotide_index=self.first_nucleotide_index+self.number_coding_nucleotides
        self.coding_nucleotide_index=numpy.array([i for i in range(self.first_nucleotide_index,self.last_nucleotide_index,1)])
        self.coding_nucleotides_position=numpy.arange(1,self.number_coding_nucleotides+1,1)

        self.number_promoter_nucleotides=len(self.promoter_nucleotides_string)
        self.promoter_nucleotides_position=numpy.arange(-1*(self.number_promoter_nucleotides),0,1)

        print(self.gene_name,self.reverse) #PWF

        if self.reverse:
            self.promoter_nucleotide_index=numpy.array([i for i in range (self.last_nucleotide_index,self.last_nucleotide_index+self.number_promoter_nucleotides,1)])
        else:
            self.promoter_nucleotide_index=numpy.array([i for i in range (self.first_nucleotide_index-self.number_promoter_nucleotides,self.first_nucleotide_index,1)])

        if self.element_type in ['GENE','LOCUS']:

            self.first_amino_acid_position=first_amino_acid_position

            # insist that the number of nucleotides must be divisible by three so can be decoded into amino acids
            assert (self.number_coding_nucleotides%3)==0, "string of nucleotides must be exactly divisible by three"

            # how many amino acids does that therefore correspond to?
            self.number_amino_acids=int(self.number_coding_nucleotides/3)

            self.last_amino_acid_position=self.first_amino_acid_position+self.number_amino_acids

            self._setup_conversion_dicts()

            # store the positions of the amino acids, as setup by first_amino_acid_position
            self.amino_acid_position=numpy.array([i for i in range(self.first_amino_acid_position,self.last_amino_acid_position,1)])

        self.mutation_ref=[]
        self.mutation_pos=[]
        self.mutation_new=[]

        # lastly, recompute the sequence
        self._recompute_sequence()

    def _recompute_sequence(self):

        # create numpy arrays of the bases
        self.coding_nucleotides=numpy.array([i for i in self.coding_nucleotides_string])
        self.promoter_nucleotides=numpy.array([i for i in self.promoter_nucleotides_string])

        # first find out the reverse complement if strand==-1 in the GenBank file
        if self.reverse:
            self.coding_string=self._reverse_complement(self.coding_nucleotides_string)
            self.coding_sequence=numpy.array([i for i in self.coding_string])

            rev_nucs=self._reverse_complement(self.promoter_nucleotides_string)
            self.promoter_sequence=numpy.array([i for i in rev_nucs])
        else:
            self.coding_string=self.coding_nucleotides_string
            self.coding_sequence=numpy.array([i for i in self.coding_string])
            self.promoter_sequence=numpy.array([i for i in self.promoter_nucleotides_string])

        # only translate the sequence if it is a protein coding element
        if self.element_type in ['GENE','LOCUS']:
            self._translate_sequence()

    def _translate_sequence(self):

        self.triplets=numpy.array([self.coding_string[i:i+3] for i in range(0,len(self.coding_string),3)])

        # now translate the triplets into amino acids using this new dictionary
        self.amino_acids=numpy.array([self.triplet_to_amino_acid[i] for i in self.triplets])

        # and build the amino acid string
        self.amino_acid_string=""
        for i in self.amino_acids:
            self.amino_acid_string+=i


    def _setup_conversion_dicts(self):

        bases = ['t', 'c', 'a', 'g', 'x', 'z']
        aminoacids = 'FFLLXZSSSSXZYY!!XZCC!WXZXXXXXXZZZZXZLLLLXZPPPPXZHHQQXZRRRRXZXXXXXXZZZZXZIIIMXZTTTTXZNNKKXZSSRRXZXXXXXXZZZZXZVVVVXZAAAAXZDDEEXZGGGGXZXXXXXXZZZZXZXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXZZZZXZZZZZXZZZZZXZZZZZXZXXXXXXZZZZXZ'
        self.codons = numpy.array([a+b+c for a in bases for b in bases for c in bases])
        self.triplet_to_amino_acid = dict(zip(self.codons, aminoacids))
        self.amino_acids_of_codons=numpy.array([self.triplet_to_amino_acid[i] for i in self.codons])

    def __repr__(self):

        output=self.gene_name+" gene\n"
        output+=self.element_type.upper()+"\n"
        output+="%i nucleotides\n" % self.number_coding_nucleotides
        output+=self.coding_nucleotides_string[0:13]+"..."+self.coding_nucleotides_string[-16:]+"\n"
        output+="%i to %i\n" % (self.coding_nucleotide_index[0],self.coding_nucleotide_index[-1])
        if self.element_type in ['GENE','LOCUS']:
            output+="%i nucleotides promoter\n" % (self.number_promoter_nucleotides)
            output+="%i amino acids\n" % self.number_amino_acids
            output+=self.amino_acid_string[0:4]+"..."+self.amino_acid_string[-5:]+"\n"
            output+="%i to %i\n" % (self.amino_acid_position[0],self.amino_acid_position[-1])
        return(output)

    def store_indel(self, mutation, ref=None, alt=None):

        cols=mutation.split("_")

        assert len(cols)==3, "indel passed not in form 1300_ins_3"

        self.indels[mutation]={ 'REF':ref.lower(), 'ALT':alt.lower() }

    def mutate_base(self,position=None,original_base=None,new_base=None):

        assert original_base.lower() in ['a','c','t','g'], "not been given a nucleotide!"

        assert position>-101, "position not an integer >= -100 (for promoter)"

        assert new_base.lower() in ['a','c','t','g','x','z'], "not been given a nucleotide!"

        assert original_base!=new_base, "not a mutation!"

        # is the position in the coding sequence or the promoter?
        if (position in self.coding_nucleotide_index):

            # make a Boolean mask identifying the position to mutate
            location=self.coding_nucleotide_index==position

            # more paranoia
            assert numpy.sum(location)==1, "WARNING: trying to mutate "+self.gene_name+" at position "+str(position)+" from "+original_base+" to "+new_base+" and the mask has "+str(numpy.sum(location))+" locations"

            # check that the passed reference base matches what we think it should be
            if self.coding_nucleotides[location][0]!=original_base.lower():
                logging.warning(self.gene_name+","+str(position)+",supplied wildtype base "+original_base.lower()+" does not match coding of "+self.coding_nucleotides[location][0]+" in the GenBank file! (Genbank version mismatch?)" )

            # mutate to the new base
            self.coding_nucleotides[location]=new_base.lower()

            # recreate the string
            self.coding_nucleotides_string=""
            for i in self.coding_nucleotides:
                self.coding_nucleotides_string+=i

        elif position in self.promoter_nucleotide_index:

            # make a Boolean mask identifying the position to mutate
            location=self.promoter_nucleotide_index==position

            # more paranoia
            assert numpy.sum(location)==1, "WARNING: trying to mutate "+self.gene_name+" at position "+str(position)+" from "+original_base+" to "+new_base+" and the mask has "+str(numpy.sum(location))+" locations"

            # check that it is actually a mutation!
            assert self.promoter_nucleotides[location]!=new_base.lower(),  "new base is also "+new_base

            # mutate to the new base
            self.promoter_nucleotides[location]=new_base.lower()

            # recreate the string
            self.promoter_nucleotides_string=""
            for i in self.promoter_nucleotides:
                self.promoter_nucleotides_string+=i

        else:

            logging.warning(self.gene_name+","+str(position)+",position not in CDS or promoter of gene")

        self._recompute_sequence()

    def _reverse_complement(self,nucleotides_string):

        complementary_bases = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'x':'x', 'z':'z'}

        reversed_nucleotides_string=nucleotides_string[::-1]

        bases=list(reversed_nucleotides_string)

        reverse_complement=""
        for i in [complementary_bases.get(base,base) for base in bases]:
            reverse_complement+=i

        return(reverse_complement)

    def return_location(self,position):

        if position in self.coding_nucleotide_index:
            return("CDS",((position-self.coding_nucleotide_index[0])+1))
        elif position in self.promoter_nucleotide_index:
            return("PROMOTER",position-self.coding_nucleotide_index[0])
        else:
            logging.warning(self.gene_name+","+str(position)+",supplied position not in gene ")
            return(None,None)


    def identify_mutations(self,reference):

        # let's check we've been given a reference gene with the same name and same number of nucleotides
        assert self.gene_name==reference.gene_name
        assert self.number_coding_nucleotides==reference.number_coding_nucleotides
        assert self.number_promoter_nucleotides==reference.number_promoter_nucleotides

        self.mutations={}

        promoter_mutation_mask=self.promoter_sequence!=reference.promoter_sequence

        # self.mutation_ref=reference.promoter_nucleotides[promoter_mutation_mask]
        # self.mutation_new=self.promoter_nucleotides[promoter_mutation_mask]
        # self.mutation_pos=self.promoter_nucleotides_position[promoter_mutation_mask]

# it had better not be this!
        self.mutation_ref=reference.promoter_sequence[promoter_mutation_mask]
        self.mutation_new=self.promoter_sequence[promoter_mutation_mask]
        self.mutation_pos=self.promoter_nucleotides_position[promoter_mutation_mask]

        if self.element_type in ['GENE','LOCUS']:

            # make Boolean array identifying where there are differences between the nucleotide triplets/amino acids
            coding_mutation_mask=self.triplets!=reference.triplets

            self.mutation_ref=numpy.append(self.mutation_ref,reference.triplets[coding_mutation_mask])
            self.mutation_new=numpy.append(self.mutation_new,self.triplets[coding_mutation_mask])
            self.mutation_pos=numpy.append(self.mutation_pos,self.amino_acid_position[coding_mutation_mask])

        elif self.element_type=="RNA":

            coding_mutation_mask=self.coding_sequence!=reference.coding_sequence

            self.mutation_ref=numpy.append(self.mutation_ref,reference.coding_nucleotides[coding_mutation_mask])
            self.mutation_new=numpy.append(self.mutation_new,self.coding_nucleotides[coding_mutation_mask])
            self.mutation_pos=numpy.append(self.mutation_pos,self.coding_nucleotides_position[coding_mutation_mask])

        if self.gene_name=="Rv1258c":
            print("promoter_nucleotides")
            print(reference.promoter_nucleotides)
            print(self.promoter_nucleotides)
            print(reference.promoter_nucleotides[promoter_mutation_mask])
            print(self.promoter_nucleotides[promoter_mutation_mask])
            print("promoter_sequence")
            print(reference.promoter_sequence)
            print(self.promoter_sequence)
            print(reference.promoter_sequence[promoter_mutation_mask])
            print(self.promoter_sequence[promoter_mutation_mask])
            print("final")
            print(self.mutation_ref)
            print(self.mutation_new)
            print(self.mutation_pos)

        for i in self.indels:

            cols=i.split("_")

            assert len(cols)==3, "malformed indel mutation"

            if int(cols[0])<0:
                promoter=True
                cds=False
            else:
                cds=True
                promoter=False

            if int(cols[2])>0:
                insertion=True
                deletion=False
            else:
                insertion=False
                deletion=True

            if cols[2]=="*":
                number_nucleotide_changes=None
            else:
                number_nucleotide_changes=int(cols[2])

            self.mutations[i]={ "ELEMENT_TYPE": self.element_type,\
                                "VARIANT_TYPE": "INDEL",
                                "POSITION": int(cols[0]),\
                                "PROMOTER":promoter,
                                "CDS":cds,\
                                "SYNONYMOUS":False,\
                                "NONSYNONYMOUS":False,\
                                "INSERTION":insertion,\
                                "DELETION":deletion,\
                                "REF":self.indels[i]['REF'],\
                                "ALT":self.indels[i]['ALT'],\
                                "NUMBER_NUCLEOTIDE_CHANGES":number_nucleotide_changes}


        for (before,position,after) in zip(self.mutation_ref,self.mutation_pos,self.mutation_new):

            if (self.element_type=='RNA') or (self.element_type in ['GENE','LOCUS'] and position<0):

                final_before=before
                final_after=after
                # before=None
                # after=None
                number_nucleotide_changes=1

                assert before!=after

            else:

                # count the number of nucleotides different between the two triplets
                number_nucleotide_changes=sum(1 for a, b in zip(before, after) if a != b)

                # convert the triplets to amino acids
                final_before=self.triplet_to_amino_acid[before]
                final_after=self.triplet_to_amino_acid[after]

            if position<0:
                cds=False
                promoter=True
            else:
                cds=True
                promoter=False

            mutation_name=final_before+str(position)+final_after

            if final_before==final_after:
                synonymous=True
                nonsynonymous=False
            else:
                synonymous=False
                nonsynonymous=True

            self.mutations[mutation_name]={ "ELEMENT_TYPE": self.element_type,\
                                "VARIANT_TYPE": "SNP",
                                "POSITION": position,\
                                "PROMOTER":promoter,
                                "CDS":cds,\
                                "SYNONYMOUS":synonymous,\
                                "NONSYNONYMOUS":nonsynonymous,\
                                "INSERTION":False,\
                                "DELETION":False,\
                                "REF":before,\
                                "ALT":after,\
                                "NUMBER_NUCLEOTIDE_CHANGES":number_nucleotide_changes}


    def print_mutations(self):

        line=""
        for (i,j,k) in zip(self.mutation_ref,self.mutation_pos,self.mutation_new):
            line+=("%s%s%s " % (i,j,k))
        return line
