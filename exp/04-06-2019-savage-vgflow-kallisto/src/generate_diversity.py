# Generating simulated pathogen diversity
# Cedric Chauve, cedric.chauve@sfu.ca

import sys
import random
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

## ----------------------------------------------------------------------
# Reading Strain types FASTA sequences

# Output:
# dictionary of sequences indexed by ST
def read_ST_sequences(file_name):
    data = open(file_name,'r').readlines()
    STSequencesDict = {}
    for l in data:
        l1 = l.rstrip().split('\t')
        if l1[0][0] == '>':
            st = l1[0].replace('>','')
        else:
            STSequencesDict[st] = l.rstrip()
    return(STSequencesDict)

## ----------------------------------------------------------------------
# Generating random diversity configuration

# Select nb_strains random STs, each with a copy number at most max_copy_number
def simulate_diversity(ST_sequences_dict, nb_strains, max_copy_number):
    STList   = list(ST_sequences_dict.keys())
    STChosen      = random.sample(STList,nb_strains)
    STCopyNumbers = [random.randint(1,max_copy_number) for st in range(nb_strains)]
    return({STChosen[i]: STCopyNumbers[i]  for i in range(nb_strains)})

# Creating a fasta file for a random diversity, accounting for copy number
def diversity_2_fasta(diversity,ST_sequences_dict,output_name):
    output = open(output_name,'w')
    # Creating fasta sequences per ST
    STList = diversity.keys() # Strain types
    for st in STList:
        cn         = diversity[st] # Copy number of the ST
        STSequence = ST_sequences_dict[st]   # Alleles of the ST
        # Writing in the FASTA file
        for c in range(cn):
            output.write('>'+st+'_'+str(c+1)+'\n'+STSequence+'\n')

## ----------------------------------------------------------------------
if __name__ == "__main__":
    ST_FILE     = sys.argv[1]
    NB_STRAINS  = int(sys.argv[2])
    MAX_COPY_NB = int(sys.argv[3])
    SEED        = int(sys.argv[4])
    OUTPUT_DIR  = sys.argv[5]
    PREFIX      = sys.argv[6]

    LogName =  OUTPUT_DIR+'/'+PREFIX+'.log'
    Log     = open(LogName,'w')
    
    random.seed(SEED)

    STSequencesDict = read_ST_sequences(ST_FILE)
    Diversity       = simulate_diversity(STSequencesDict,NB_STRAINS,MAX_COPY_NB)
    Log.write('# Random configuration with seed '+str(SEED)+'\n')
    for st in Diversity.keys():
        Log.write(st+'\t'+str(Diversity[st])+'\n')
    DiversityFasta = OUTPUT_DIR+'/'+PREFIX+'.fasta'
    diversity_2_fasta(Diversity,STSequencesDict,DiversityFasta)
    Log.close()
