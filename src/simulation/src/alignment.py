import os
import subprocess

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def align(strains):
  result = None
  seqs   = (SeqRecord(s.concatenated(), id=s.id) for s in strains)
  
  with subprocess.Popen(
    [os.environ["MUSCLE"]],
    stdin    = subprocess.PIPE,
    stdout   = subprocess.PIPE,
    stderr   = subprocess.DEVNULL,
    encoding = "ascii"
  ) as proc:
    SeqIO.write(seqs, proc.stdin, "fasta")
    proc.stdin.close()
    result = {s.id : s.seq for s in SeqIO.parse(proc.stdout, "fasta")}
    
  return result

# r = align([Seq("ATCGATCGATCG"), Seq("TCGATATCG")])
# for s in r:
#   print(str(s))
  
