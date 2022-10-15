#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HOMOGENOUES_MODELS

Given aligned reads from a disomy and a monosomy, the likelihood of observed
reads under the unmatched and matched hypothesis is calculated. Here we assume
a non-admixed ancestral makeup.

The matched hypothesis means that the disomy and the monosomy have matched
haplotypes. In contrast, under the unmatched hypothesis the disomy and the
monosomy do not share common haplotypes.

Daniel Ariad (daniel@ariad.org)
Dec 21, 2020
"""

import pickle, os, sys, bz2, collections, gzip, platform

from functools import reduce
from operator import and_, itemgetter
from itertools import combinations

leg_tuple = collections.namedtuple('leg_tuple', ('chr_id', 'pos', 'ref', 'alt')) #Encodes the rows of the legend table
obs_tuple = collections.namedtuple('obs_tuple', ('pos', 'read_id', 'base')) #Encodes the rows of the observations table

### Getting a function to count non-zero bits in positive integer.
class popcount_lk:
    """ Creates an instance for calculating the population count of
        bitstring, based on a lookup table of 8 bits. """

    def __init__(self):
        """ Creates a large lookup table of the Hamming weight of every 8 bit integer. """
        self.lookup_table = bytes.maketrans(bytes(range(1<<8)),bytes((bin(i).count('1') for i in range(1<<8))))
        self.byteorder = sys.byteorder

    def __call__(self,x):
        """ Breaks x, which is a python integer type, into chuncks of 8 bits.
        Calls the lookup table to get the population count of each chunck and returns
        the aggregated population count. """

        return sum(x.to_bytes((x.bit_length()>>3)+1,self.byteorder).translate(self.lookup_table))

class homogeneous:
    """ Based on the statisitcal models for non-admixed and a reference
    panel, it allows to calculate the likelihoods of observed alleles under the
    MATCHED and UNMATCHED statistical models. """


    def __init__(self, obs_tab, leg_tab, hap_tab_per_group, number_of_haplotypes_per_group):
        """ Initialize the attributes of the class. """


        print('Initializing the algorithm for non-admixtures:')

        if len(hap_tab_per_group)!=1:
            raise Exception('Error: The non-admixed models can not be applied when the HAP file includes multiple superpopulations/group2.')

        group2 = [*hap_tab_per_group.keys()].pop()
        hap_tab = hap_tab_per_group[group2]
        self.number_of_haplotypes = number_of_haplotypes_per_group[group2]

        if len(leg_tab)!=len(hap_tab):
            raise Exception('Error: the number of SNPs in the LEGEND file differ from the number of SNPs in the HAP file.')


        self.fraction_of_matches = {}
        self.hap_dict, self.fraction_of_matches[group2] = self.build_hap_dict(obs_tab, leg_tab, hap_tab, self.number_of_haplotypes)


        print('%.2f%% of the observed alleles matched the reference panel.' % (100*self.fraction_of_matches[group2]))

    def build_hap_dict(self, obs_tab, leg_tab, hap_tab, number_of_haplotypes):
        """ Returns a dictionary that lists SNP alleles and gives their
        relevent row from haplotypes table. The row is stored as bits, where
        True means that the haplotype contains the allele. We denote the
        returned dictionary as the reference panel. """

        hap_dict = dict()
        mismatches = 0
        combined = {pos: (ref,alt,hap) for (chr_id,pos,ref,alt),hap in zip(leg_tab, hap_tab)}
        missing = 3*(None,)

        b = (1 << number_of_haplotypes) - 1 #### equivalent to int('1'*number_of_haplotypes,2)

        for (pos, read_id, base) in obs_tab:
            ref, alt, hap = combined.get(pos, missing)
            if base==alt:
                hap_dict[(pos,base)] = hap
            elif base==ref:
                hap_dict[(pos,base)] = hap ^ b ### ^b flips all bits of the binary number, hap_tab[ind] using bitwise xor operator.
            else:
                mismatches += 1

        fraction_of_matches = 1-mismatches/len(obs_tab)

        return hap_dict, fraction_of_matches

    def build_intrenal_hap_dict(self, alleles):
        """ This function allows treatment of alleles and haplotypes on an
        equal footing. This is done in three steps: (1) All the alleles and
        haplotypes are enumerated. (2) For each given haplotype, tuples in the
        reference panel, associated with the haplotype's alleles, are
        intersected. (3) A dictionary that lists all the alleles and haplotypes
        by their index is returned. The dictionary gives for each allele and
        haplotype their associated tuple and intersected tuple, respectively. """

        internal = {}
        for i, haplotype in enumerate(alleles):
            if type(haplotype[0])==tuple: #Checks if X is a tuple/list of alleles.
                n = len(haplotype)
                if n==1:
                    internal[1 << i] = self.hap_dict[haplotype[0]]
                elif n==2:
                    internal[1 << i] = self.hap_dict[haplotype[0]] & self.hap_dict[haplotype[1]]
                else:
                    internal[1 << i] = reduce(and_,itemgetter(*haplotype)(self.hap_dict))

            elif type(haplotype[0])==int: #Checks if X is a single allele.
                internal[1 << i] = self.hap_dict[haplotype]
            else:
                raise Exception('error: joint_frequencies only accepts alleles and tuple/list of alleles.')

        return internal

    def joint_frequencies_combo(self, alleles, normalize):
        """ Based on the reference panel, it calculates joint frequencies of
            observed alleles. The function arguments are alleles, that is,
            tuples of position and base, e.g., (100,'T'), (123, 'A') and
            (386, 'C'). Each allele is enumerated according to the order it
            was received by the function. The function returns a dictionary that
            lists all the possible subgroups of the given alleles. Each key in
            the dictionary is a tuple of intergers that are lexicographically
            sorted. Moreover, each integer within the keys corresponds to an
            enumerated allele. For each subgroup of alleles the dictionary
            gives the joint frequencies in a given population. The function
            arguments can also include haplotypes, that is, tuples of alleles;
            Haplotypes are treated in the same manner as alleles. """

        hap = self.build_intrenal_hap_dict(alleles)

        result = {c: popcount(A) for c,A in hap.items() }

        for C in combinations(hap, 2):
            result[C[0]|C[1]] = popcount(hap[C[0]]&hap[C[1]])

        for C in combinations(hap, 3):
            result[C[0]|C[1]|C[2]] = popcount(hap[C[0]]&hap[C[1]]&hap[C[2]])

        for r in range(4,len(alleles)):
            for C in combinations(hap, r):
                result[sum(C)] = popcount(reduce(and_,itemgetter(*C)(hap)))

        if len(alleles)>=4:
            result[sum(hap.keys())] = popcount(reduce(and_,hap.values()))

        if normalize:
            result = {k: v/self.number_of_haplotypes for k,v in result.items()}

        return result

    def joint_frequencies_redundant(self, alleles, normalize):
        """ Performs the same task as the function joint_frequencies_combo, but
            sacrifice the efficency in favor of simplicity. """

        internal = {1 << i: reduce(and_,itemgetter(*haplotype)(self.hap_dict))
                        if type(haplotype[0])==tuple else self.hap_dict[haplotype]
                                for i, haplotype in enumerate(alleles)}

        result = {c: popcount(A) for c,A in internal.items() }

        result |= {sum(C): popcount(reduce(and_,itemgetter(*C)(internal)))
                      for r in range(1,len(alleles))
                              for C in combinations(internal, r+1)}

        if normalize:
            result = {k: v/self.number_of_haplotypes
                          for k,v in result.items()}

        return result

    def likelihoods(self, alleles):
        """ Calculates the likelihood to observe a set with alleles and
        haplotypes under MATCHED and UNMATCHED scenarios. """

        F = self.joint_frequencies_combo(alleles, normalize=False)
        N = self.number_of_haplotypes #Divide values by N to normalize the joint frequencies.

        F[0] = N  ### The unnormalized frequency of a null set of alleles is N.

        L = len(alleles) - 1  ### Number of alleles observed in the disomy
        MONOSOMY = 1 << L ### The last read in the tuple is from the monosomy.
        b = (1 << L) - 1 ### A "bit array" of length L that is filled with 1-bits, i.e., 0b11111....1.

        ### UNMATCHED ###
        A =  1 / ((1 << (L-1)) * N**2)
        DISOMY = A * sum(F[B0] * F[B0 ^ b] for B0 in range(1 << (L-1)))
        UNMATCHED = F[MONOSOMY] * DISOMY / N

        ### MATCHED ###
        A =  1 / ((1 << L) * N**2)
        MATCHED = A * sum(F[B0] * F[(B1 := B0 ^ b)|MONOSOMY] + F[B1] * F[B0|MONOSOMY]
                          for B0 in range(1 << (L-1)))

        return UNMATCHED, MATCHED


    def likelihoods2(self, alleles):
        """ Calculates the likelihood to observe two alleles/haplotypes
        under MATCHED and UNMATCHED scenarios. """

        F = self.joint_frequencies_combo(alleles, normalize=True)
        a, b, ab = F[1], F[2], F[3]

        UNMATCHED = a*b
        MATCHED = (ab+a*b)/2

        return UNMATCHED, MATCHED

    def likelihoods3(self, alleles):
        """ Calculates the likelihood to observe three alleles/haplotypes
        under MATCHED and UNMATCHED scenarios. """

        F = self.joint_frequencies_combo(alleles, normalize=True)
        a, b, c = F[1], F[2], F[4]
        ab, ac, bc =  F[3], F[5], F[6]
        abc = F[7]

        UNMATCHED = c*(ab+a*b)/2
        MATCHED = (ab*c+ac*b+abc+a*bc)/4
        return UNMATCHED, MATCHED

    def likelihoods4(self, alleles):
        """ Calculates the likelihood to observe four alleles/haplotypes
        under MATCHED and UNMATCHED scenarios. """

        F = self.joint_frequencies_combo(alleles, normalize=True)
        a, b, c, d = F[1], F[2], F[4], F[8]
        ab, ac, ad, bc, bd, cd = F[3], F[5], F[9], F[6], F[10], F[12]
        abc, abd, acd, bcd = F[7], F[11], F[13], F[14]
        abcd = F[15]

        UNMATCHED = d*(abc+ab*c+ac*b+bc*a)/4 #The likelihood of diploidy.
        MATCHED = (abc*d+ab*cd+ac*bd+bc*ad+abcd+abd*c+acd*b+bcd*a)/8
        return UNMATCHED, MATCHED

    def likelihoods5(self, alleles):
        """ Calculates the likelihood to observe five alleles/haplotypes
        under MATCHED and UNMATCHED scenarios. """

        F = self.joint_frequencies_combo(alleles, normalize=True)
        a, b, c, d, e = F[1], F[2], F[4], F[8], F[16]
        ab, ac, ad, ae, bc, bd, be, cd, ce, de = F[3], F[5], F[9], F[17], F[6], F[10], F[18], F[12], F[20], F[24]
        abc, abd, abe, acd, ace, ade, bcd, bce, bde, cde = F[7], F[11], F[19], F[13], F[21], F[25], F[14], F[22], F[26], F[28]
        abcd, abce, abde, acde, bcde = F[15], F[23], F[27], F[29], F[30]
        abcde = F[31]

        UNMATCHED = e*(abcd+abc*d+bcd*a+acd*b+abd*c+ab*cd+ad*bc+ac*bd)/8
        MATCHED = (abcd*e+abc*de+bcd*ae+acd*be+abd*ce+ab*cde+ad*bce+ac*bde+abcde+abce*d+bcde*a+acde*b+abde*c+abe*cd+ade*bc+ace*bd)/16 #The likelihood of diploidy. #V

        return UNMATCHED, MATCHED

    def get_likelihoods(self, alleles):
        """ Uses the optimal function to calculate the likelihoods.
        In general, self.likelihoods can get more than five alleles but the
        dedicated functions are optimized to a certain number of alleles. """

        l = len(alleles)
        if l==2:
            result = self.likelihoods2(alleles)
        elif l==3:
            result = self.likelihoods3(alleles)
        elif l==4:
            result = self.likelihoods4(alleles)
        elif l==5:
            result = self.likelihoods5(alleles)
        else:
            result = self.likelihoods(alleles)
        return result

def wrapper_of_homogenoues_for_debugging(obs_filename,leg_filename,hap_filename,group2):
    """ Wrapper function of the class 'homogeneous'. It receives an observations
    file, legend file, haplotypes file, samples file and a file with the
    statistical models. Based on the given data it creates and returns an
    instance of the class. """

    if not os.path.isfile(obs_filename): raise Exception('Error: OBS file does not exist.')
    if not os.path.isfile(leg_filename): raise Exception('Error: LEGEND file does not exist.')
    if not os.path.isfile(hap_filename): raise Exception('Error: HAP file does not exist.')

    load = lambda filename: {'bz2': bz2.open, 'gz': gzip.open}.get(filename.rsplit('.',1)[1], open)  #Adjusts the opening method according to the file extension.

    open_hap = load(hap_filename)
    with open_hap(hap_filename,'rb') as hap_in:
        hap_tab_per_group = pickle.load(hap_in)
        number_of_haplotypes_per_group = pickle.load(hap_in)

    open_leg = load(leg_filename)
    with open_leg(leg_filename,'rb') as leg_in:
        leg_tab = pickle.load(leg_in)

    open_obs = load(obs_filename)
    with open_obs(obs_filename, 'rb') as obs_in:
        obs_tab = pickle.load(obs_in)
        #info = pickle.load(f)

    return homogeneous(obs_tab, leg_tab, {group2: hap_tab_per_group[group2]}, {group2: number_of_haplotypes_per_group[group2]})

if platform.python_implementation()=='PyPy':
    popcount = popcount_lk()
else:
    try:
        from gmpy2 import popcount
    except Exception as error_msg:
        print(error_msg)
        popcount = popcount_lk()

if __name__ != "__main__":
    print('The module HOMOGENOUES_MODELS was imported.')
else:
    print('The module HOMOGENOUES_MODELS was invoked directly')
    sys.exit(0)

###############################   END OF FILE   ###############################



"""
if __name__ != "__main__":
    print("The module HOMOGENOUES_MODELS was imported.")
else:
    print("The module HOMOGENOUES_MODELS was invoked directly")
    #sys.exit(0)
    import time, random
    t0 = time.time()
    obs_filename = 'test/test.obs.p'
    hap_filename = 'test/test.hap.p'
    leg_filename = 'test/test.leg.p'

    A = wrapper_of_homogenoues_for_debugging(obs_filename,leg_filename,hap_filename,'group1')

    alleles = tuple(A.hap_dict.keys())

    frequencies = lambda x: {bin(a)[2:]:b for a,b in A.joint_frequencies_combo(x,normalize=True).items()}
    frequencies_redundant = lambda x: {bin(a)[2:]:b for a,b in A.joint_frequencies_redundant(x,normalize=True).items()}

    likelihoods = A.likelihoods
    likelihoods2 = A.likelihoods2
    likelihoods3 = A.likelihoods3
    likelihoods4 = A.likelihoods4
    likelihoods5 = A.likelihoods5


    random.seed(a=2022, version=2)
    x = random.randrange(len(alleles)-16) #123
    haplotypes = (alleles[x:x+4],alleles[x+4:x+8],alleles[x+8:x+12],alleles[x+12:x+16])

    print('-----joint_frequencies_combo-----')
    print(frequencies(alleles[x+0:x+1]))
    print(frequencies_redundant(alleles[x+0:x+1]))
    print(frequencies(alleles[x:x+4]))
    print(frequencies_redundant(alleles[x:x+4]))
    print('-----likelihoods4-haplotypes-----')
    print(haplotypes)
    print(frequencies(haplotypes))
    print(frequencies_redundant(haplotypes))
    print(likelihoods(haplotypes))
    print(likelihoods4(haplotypes))
    print('-----likelihoods2-----')
    print(alleles[x:x+2])
    print(frequencies(alleles[x:x+2]))
    print(frequencies_redundant(alleles[x:x+2]))
    print(likelihoods(alleles[x:x+2]))
    print(likelihoods2(alleles[x:x+2]))
    print('-----likelihoods3-----')
    print(alleles[x:x+3])
    print(frequencies(alleles[x:x+3]))
    print(frequencies_redundant(alleles[x:x+3]))
    print(likelihoods(alleles[x:x+3]))
    print(likelihoods3(alleles[x:x+3]))
    print('-----likelihoods4-----')
    print(alleles[x:x+4])
    print(frequencies(alleles[x:x+4]))
    print(frequencies_redundant(alleles[x:x+4]))
    print(likelihoods(alleles[x:x+4]))
    print(likelihoods4(alleles[x:x+4]))
    print('-----likelihoods5-----')
    print(alleles[x:x+5])
    print(frequencies(alleles[x:x+5]))
    print(frequencies_redundant(alleles[x:x+5]))
    print(likelihoods(alleles[x:x+5]))
    print(likelihoods5(alleles[x:x+5]))

    t1 = time.time()

    print('Done in %.3f sec.' % ((t1-t0)))

"""


"""
T = []
for p in range(1,32):
    x = ''.join([chr(ord('a')+i) for i,j in enumerate(bin(p)[:1:-1]) if j=='1'])
    T.append((x,p))
T2 = sorted(T,key=lambda x: (len(x[0]),x[0]) )
for x,p in T2:
    r = f'{x:s} = F[{p:d}]'
    print(r)
"""
