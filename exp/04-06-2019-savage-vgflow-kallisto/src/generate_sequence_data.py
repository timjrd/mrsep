# Generating simulated sequencing data of pathogen diversity
# Cedric Chauve, cedric.chauve@sfu.ca

import os
import sys
import random

# Generating simulated sequence data

ART = '../../../bin/art_bin_MountRainier/art_illumina'

def artCommand(sample_fasta,sample_name,coverage,hss,read_lg,insert_lg):
    return(ART+' -ss '+str(hss)+' -sam -i '+sample_fasta+' -p -l '+str(read_lg)+' -f '+str(coverage)+' -m '+str(insert_lg)+' -s 10 -o '+sample_name+'_')

## ----------------------------------------------------------------------
if __name__ == "__main__":
    DIV_FASTA   = sys.argv[1]
    DEPTH_COV   = int(sys.argv[2])
    READ_LG     = int(sys.argv[3])
    INSERT_LG   = int(sys.argv[4])
    SEQ_SYS     = sys.argv[5]
    OUTPUT_DIR  = sys.argv[6]
    PREFIX      = sys.argv[7]

    ArtPrefix = OUTPUT_DIR+'/'+PREFIX+'_'+SEQ_SYS
    os.system(artCommand(DIV_FASTA,ArtPrefix,DEPTH_COV,SEQ_SYS,READ_LG,INSERT_LG)+' > '+ArtPrefix+'.log')
    
