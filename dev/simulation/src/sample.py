import sys
import random

from Bio     import SeqIO
from Bio.Seq import Seq

def randomSeq(length):
  return Seq("".join([random.choice("ATCG") for _ in range(length)]))

if len(sys.argv) != 7:
  print("usage: " + sys.argv[0] + " input.fasta nb-alleles min-copy-nb max-copy-nb ends-length seed")
  sys.exit(1)

inputFile  = sys.argv[1]
nbAlleles  = int(sys.argv[2])
minCopyNb  = int(sys.argv[3])
maxCopyNb  = int(sys.argv[4])
endsLength = int(sys.argv[5])
seed       = int(sys.argv[6])

random.seed(seed)

alleles = list(SeqIO.parse(open(inputFile), "fasta"))
random.shuffle(alleles)
chosen = alleles[0:nbAlleles]

result = []
for allele in chosen:
  n = random.randint(minCopyNb, maxCopyNb)
  for i in range(n):
    copy = allele.upper()
    copy.description = ""
    copy.id = str(i) + "." + copy.id
    copy.seq = randomSeq(endsLength) + copy.seq + randomSeq(endsLength)
    result.append(copy)

random.shuffle(result)
SeqIO.write(result, sys.stdout, "fasta")
