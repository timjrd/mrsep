# Generating simulated pathogen diversity, ensuring close strains
# Cedric Chauve, cedric.chauve@sfu.ca

import sys
import random
import operator
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
# Reading Strain types distance

def read_ST_distances(distances_file):
    STDist    = {}
    distances_data = open(distances_file,'r').readlines()
    for l in distances_data:
        l1 = l.rstrip().split()
        ST = l1[0]
        STDist[ST] = {}
        for i in range(1,len(l1)):
            l2 = l1[i].split(':')
            STDist[ST][l2[0]] = int(l2[1])
            STDist[l2[0]][ST] = int(l2[1])
    return(STDist)


## ----------------------------------------------------------------------
# Generating random diversity configuration

# Select nb_strains random STs, each with a copy number at most max_copy_number
def simulate_diversity(ST_dist_dict, nb_strains, max_copy_number):
    STList    = list(ST_dist_dict.keys())
    STChosen  = random.sample(STList,1)
    STCurrent = STChosen[0]
    for i in range(1,nb_strains):
        STSortedByDist = sorted(ST_dist_dict[STCurrent].items(), key=operator.itemgetter(1), reverse=False)
        # Next strain type: the closest not already in the list
        j = 0
        STClosest = STSortedByDist[j][0]
        while STClosest in STChosen:
            j += 1
            STClosest = STSortedByDist[j][0]
        STChosen += [STClosest]
        STCurrent = STClosest
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
    ST_DIST_FILE = sys.argv[2]
    NB_STRAINS  = int(sys.argv[3])
    MAX_COPY_NB = int(sys.argv[4])
    SEED        = int(sys.argv[5])
    OUTPUT_DIR  = sys.argv[6]
    PREFIX      = sys.argv[7]

    LogName =  OUTPUT_DIR+'/'+PREFIX+'.log'
    Log     = open(LogName,'w')
    
    random.seed(SEED)

    STSequencesDict = read_ST_sequences(ST_FILE)
    STDistancesDict = read_ST_distances(ST_DIST_FILE)
    Diversity       = simulate_diversity(STDistancesDict,NB_STRAINS,MAX_COPY_NB)
    Log.write('# Random configuration with seed '+str(SEED)+'\n')
    for st in Diversity.keys():
        Log.write(st+'\t'+str(Diversity[st])+'\n')
    DiversityFasta = OUTPUT_DIR+'/'+PREFIX+'.fasta'
    diversity_2_fasta(Diversity,STSequencesDict,DiversityFasta)
    Log.close()
