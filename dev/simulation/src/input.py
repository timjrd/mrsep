import sys
from collections import defaultdict

from Bio     import SeqIO
from Bio.Seq import Seq

def stat(k,v):
  sys.stderr.write(k + " " + str(v) + "\n")

def parseAln(takeRef,stream):
  result = []
  ref    = None
  offset = None
  rev    = None
  skip   = None
  for line in stream:
    if line[0] == '>':
      words  = line[1:].split()
      ref    = "".join(words[0].split('.')[1:])
      offset = int(words[2])
      rev    = words[3] == "-"
      skip   = not takeRef
    elif skip == True:
      skip   = False
    elif skip == False:
      read   = Seq(line.strip())
      read   = read.reverse_complement() if rev else read
      result.append((ref,offset,rev,read))
      ref    = None
      offset = None
      rev    = None
      skip   = None
  return result

def span(endsLength, refs, (ref,offset,rev,read)):
  offset = offset-endsLength
  edge   = offset < 0 and abs(offset) < len(read)
  read   = (read[:offset] if rev else read[abs(offset):]) if edge else read
  offset = 0 if edge else offset
  if offset < 0:
    return []
  else:
    result   = []
    ref      = filter(lambda (_,x): x != '-', enumerate(refs[ref].seq))    
    offset   = len(ref)-len(read)-offset if rev else offset
    overflow = len(ref)-len(read)-offset
    read     = read[abs(offset):] if offset < 0 else (read[:overflow] if overflow < 0 else read)
    offset   = max(0,offset)
    for i,x in enumerate(read):
      r = ref[offset+i]
      result.append((x, r[0], r[1]))
    return result

if len(sys.argv) != 5:
  print("usage: " + sys.argv[0] + " aligned-input.fasta art-output.aln ends-length take-ref")
  sys.exit(1)

refsFile      = sys.argv[1]
artOutputFile = sys.argv[2]
endsLength    = int(sys.argv[3])
takeRef       = sys.argv[4] == "true"
      
refs  = SeqIO.to_dict(SeqIO.parse(open(refsFile), "fasta"))
reads = parseAln(takeRef, open(artOutputFile))

for ref in refs.values():
  print(ref.seq)
  
print("")

cov         = defaultdict(lambda: 0)
nbBases     = 0
nbReads     = 0
faultyBases = 0
faultyReads = 0
maxFaultyBasesInRead = 0
for i,read in enumerate(reads):
  xs = span(endsLength,refs,read)
  if xs:
    faultyBasesInRead = 0
    line = str(i+1) + " 1"
    for x,k,r in xs:
      line += " " + x + " " + str(k+1)
      if x == r:
        cov[(read[0],k)] += 1;
      else:
        faultyBasesInRead += 1
    print(line)
    nbBases += len(xs)
    nbReads += 1
    if faultyBasesInRead > 0:
      faultyBases += faultyBasesInRead
      faultyReads += 1
      if faultyBasesInRead > maxFaultyBasesInRead:
        maxFaultyBasesInRead = faultyBasesInRead

for ref,seq in refs.iteritems():
  lo = sys.maxint
  hi = 0
  s  = 0
  n  = 0
  cov_all = []
  for k,r in enumerate(seq):
    if r != '-':
      c = cov[(ref,k)]
      cov_all.append(c)
      if c < lo: lo = c
      if c > hi: hi = c
      s += c
      n += 1
  stat("cov_min_" + ref, lo)
  stat("cov_max_" + ref, hi)
  stat("cov_sum_" + ref, s)
  stat("cov_nb_"  + ref, n)
  for c in cov_all:
    stat("cov_all_" + ref, c)

stat("faulty_bases"             , faultyBases)
stat("nb_bases"                 , nbBases)
stat("faulty_reads"             , faultyReads)
stat("nb_reads"                 , nbReads)
stat("max_faulty_bases_per_read", maxFaultyBasesInRead)
