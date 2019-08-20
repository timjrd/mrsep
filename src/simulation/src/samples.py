import copy
import random
import mutate as m

import Strain

from Bio.SeqRecord import SeqRecord

def _duplicate(strains):
  result = []
  for s in strains:
    n = random.randint(1,5)
    for _ in range(n):
      result.append(copy.deepcopy(s))
  return result

def samples(strains):
  st = copy.deepcopy(strains)
  random.shuffle(st)
  
  result = []
  
  k = random.randint(1,4)
  result.append(("rand", _duplicate(st[:k])))
  st = st[k:]
  
  for name, pattern in m.PATTERNS:
    sample = []
    s2 = st[:2]
    st = st[2:]
    for f,s in zip(pattern,s2):
      sample += map(lambda h: s.from_history(h), f(s.concatenated()))
    
    result.append((name, _duplicate(sample)))
  
  return result

def to_json(sample):
  return {
    "name"   : sample[0],
    "strains": list(map(lambda y: y.to_json(), sample[1]))
  }

def from_json(data):
  return [
    data["name"],
    list(map(Strain.from_json, data["strains"]))
  ]

def to_fasta(sample):
  result = []

  for strain in sample[1]:
    result.append(SeqRecord(strain.interleaved(), id=strain.id))

  return result
