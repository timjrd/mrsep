# Generating simulated pathogen diversity
# Cedric Chauve, cedric.chauve@sfu.ca

import sys
import random
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

## ----------------------------------------------------------------------
# Reading MLST scheme data: strain types and alleles sequences

# Reading the list of strain types
# Format of tab-separated input file:
# ST	clpA	clpX	nifS	pepX	pyrG	recG	rplB	uvrA
# 1	1	1	1	1	1	1	1	1	
# 2	1	1	7	1	1	1	1	11	
# 3	4	1	1	1	1	6	1	7
# ...
# Output:
# dictionary of loci indexed by their order in the input file
# dictionary indexed by ST and pointing to a dictionary indexed by loci
def read_strain_types(file_name):
    data = open(file_name,'r').readlines()
    for l in data:
        l1 = l.rstrip().split('\t')
        if l1[0] == 'ST':
            STDict = {}
            lociDict = {i: l1[i] for i in range(1,len(l1))}
        else:
            STDict[l1[0]] = {lociDict[i]: lociDict[i]+'_'+l1[i] for i in range(1,len(l1))}
    return((lociDict,STDict))

# Reading the FASTA sequences of all alleles
# Input: directory to fasta file, dictionary of loci, suffix of fasta files
# Output: Dictionary of dictionary of Seq records indexed by loci names
def read_alleles_sequences(fasta_dir,fasta_suffix,loci_dict):
    lociList      = list(loci_dict.values())
    lociSequences = {locus: {} for locus in lociList}
    for locus in lociList:
        with open(fasta_dir+'/'+locus+'.'+fasta_suffix) as handle:
            for allele in SimpleFastaParser(handle):
                lociSequences[locus][allele[0]] = allele[1]
    #lociSequences = {locus: SeqIO.to_dict(SeqIO.parse(fasta_dir+'/'+locus+'.'+fasta_suffix,"fasta")) for locus in lociList}
    return(lociSequences)

# Creating a fasta file of the STs
# Output a dictionary indexed by strain types with the sequences of the STs
def create_ST_fasta(ST_dict, loci_dict, loci_sequences, spacer_lg, output_name):
    output = open(output_name,'w')
    # List of loci in the MLST scheme
    lociList = list(loci_sequences.keys())
    # Creating one random spacer between any pair ofd successive loci
    randomSpacers = {i: ''.join([random.choice('AGTC') for x in range(spacer_lg)]) for i in range(len(lociList)+1)}
    STList = ST_dict.keys() 
    for st in STList:
        STAlleles          = ST_dict[st]   # Alleles of the ST
        STAllelesSequences = {locus: loci_sequences[locus][STAlleles[locus]] for locus in lociList} # Alleles sequences
        # Creating a single FASTA sequence for the ST
        STSequence = randomSpacers[0]
        for i in range(len(lociList)):
            locus = lociList[i]
            STSequence += STAllelesSequences[locus]+randomSpacers[i+1]
        # Writing in the FASTA file
        output.write('>'+st+'\n'+STSequence+'\n')
        
## ----------------------------------------------------------------------
if __name__ == "__main__":
    ST_FILE     = sys.argv[1]
    FASTA_DIR   = sys.argv[2]
    SPACER_LG   = int(sys.argv[3])
    SEED        = int(sys.argv[4])
    OUTPUT_DIR  = sys.argv[5]

    random.seed(SEED)
    (LociDict,STDict) = read_strain_types(ST_FILE)
    LociSequences     = read_alleles_sequences(FASTA_DIR,'fasta',LociDict)
    create_ST_fasta(STDict,LociDict,LociSequences,SPACER_LG,OUTPUT_DIR+'/ST'+str(SPACER_LG)+'.fasta'
)
