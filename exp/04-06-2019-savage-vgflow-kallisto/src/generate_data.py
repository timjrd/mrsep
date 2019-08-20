# Generating simulated sequencing data of pathogen diversity
# Cedric Chauve, cedric.chauve@sfu.ca

import os
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

## ----------------------------------------------------------------------
# Generating random diversity configuration

# Select nb_strains random STs, each with a copy number at most max_copy_number
def simulate_diversity(ST_dict, loci_dict, nb_strains, max_copy_number):
    STList   = list(ST_dict.keys())
    STChosen      = random.sample(STList,nb_strains)
    STCopyNumbers = [random.randint(1,max_copy_number) for st in range(nb_strains)]
    return({STChosen[i]: STCopyNumbers[i]  for i in range(nb_strains)})

# Creating a fasta file for a random diversity, accounting for copy number
def diversity_2_fasta(diversity,ST_dict,loci_sequences,spacer_lg,output_name):
    output = open(output_name,'w')
    # List of loci in the MLST scheme
    lociList = list(loci_sequences.keys())
    # Creating one random spacer between any pair ofd successive loci
    randomSpacers = {i: ''.join([random.choice('AGTC') for x in range(spacer_lg)]) for i in range(len(lociList)+1)}
    # Creating fasta sequences per ST
    STList = diversity.keys() # Strain types
    for st in STList:
        cn                 = diversity[st] # Copy number of the ST
        STAlleles          = ST_dict[st]   # Alleles of the ST
        STAllelesSequences = {locus: loci_sequences[locus][STAlleles[locus]] for locus in lociList} # Alleles sequences
        # Creating a single FASTA sequence for the ST
        STSequence = randomSpacers[0]
        for i in range(len(lociList)):
            locus = lociList[i]
            STSequence += STAllelesSequences[locus]+randomSpacers[i+1]
        # Writing in the FASTA file
        for c in range(cn):
            output.write('>'+st+'_'+str(c+1)+'\n'+STSequence+'\n')

## ----------------------------------------------------------------------
if __name__ == "__main__":
    ST_file     = sys.argv[1]
    FASTA_DIR   = sys.argv[2]
    NB_STRAINS  = int(sys.argv[3])
    MAX_COPY_NB = int(sys.argv[4])
    SEED        = int(sys.argv[9])
    OUTPUT_DIR  = sys.argv[10]
    PREFIX      = sys.argv[11]

    LogName =  OUTPUT_DIR+'/'+PREFIX+'.log'
    Log     = open(LogName,'w')
    
    random.seed(SEED)

    (LociDict,STDict) = read_strain_types(ST_file)
    LociSequences     = read_alleles_sequences(FASTA_DIR,'fasta',LociDict)

    Diversity         = simulate_diversity(STDict,LociDict,NB_STRAINS,MAX_COPY_NB)
    Log.write('# Random configuration with seed '+str(SEED)+'\n')
    for st in Diversity.keys():
        Log.write(st+'\t'+str(Diversity[st])+'\n')
    diversityFasta = OUTPUT_DIR+'/'+PREFIX+'.fasta'
    diversity_2_fasta(Diversity,STDict,LociSequences,INSERT_LG,diversityFasta)
    Log.close()


    DEPTH_COV   = int(sys.argv[5])
    READ_LG     = int(sys.argv[6])
    INSERT_LG   = int(sys.argv[7])
    SEQ_SYS     = argv[8]

## ----------------------------------------------------------------------
# Generating simulated sequence data

ART = '../../../bin/art_bin_MountRainier/art_illumina'

def artCommand(sample_fasta,sample_name,coverage,hss,read_lg,insert_lg):
    return(ART+' -ss '+str(hss)+' -sam -i '+sample_fasta+' -p -l '+str(read_lg)+' -f '+str(coverage)+' -m '+str(insert_lg)+' -s 10 -o '+sample_name+'_'+str(hss)+'_')

    os.system(artCommand(diversityFasta,OUTPUT_DIR+'/'+PREFIX,DEPTH_COV,SEQ_SYS,READ_LG,INSERT_LG)+' >> '+LogName)
              
