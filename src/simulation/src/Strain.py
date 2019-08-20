import copy
import json
import hashlib
import random

from Bio.Seq import Seq

def _random_seq(length):
  return Seq("".join([
    random.choice("ATCG") for _ in range(length)
  ]))

def _random_seqs(n, length):
  result = []
  for _ in range(n):
    result.append(_random_seq(length))
  return result

RANDOM_SEQS_LENGTH = 150
RANDOM_SEQS = _random_seqs(100, RANDOM_SEQS_LENGTH)

def from_json(data):
  result          = Strain("",[])
  result.id       = data["id"]
  result._alleles = list(map(Seq, data["alleles"]))
  result.history  = list(map(Seq, data["history"]))
  return result

class Strain:  
  def __init__(self, id, alleles):
    self.id = id
    self._alleles = copy.deepcopy(alleles)
    self.history = []

  def _hash_id(self):
    self.id = ""
    s = json.dumps(self.to_json(), sort_keys=True)
    self.id = "hash_" + hashlib.sha256(s.encode("utf-8")).hexdigest()[:10]
    
  def to_json(self):
    return {
      "id"     : self.id,
      "alleles": list(map(str, self._alleles)),
      "history": list(map(str, self.history))
    }
    
  def from_history(self,hist):
    result = self.from_concatenated(hist[0])
    result.history = copy.deepcopy(hist)
    result._hash_id()
    return result
    
  def from_concatenated(self,seq):
    result = copy.deepcopy(self)

    #print(len(seq))
    #print(sum(map(len,result._alleles)))
    
    if len(seq) != sum(map(len,result._alleles)):
      raise IndexError
    
    for i in range(len(result._alleles)):
      l = len(result._alleles[i])
      result._alleles[i] = seq[:l]
      seq = seq[l:]
    
    result._hash_id()
    return result
  
  def interleaved(self):
    result = copy.deepcopy(RANDOM_SEQS[0])
    i = 1
    for allele in self._alleles:
      result += allele
      result += RANDOM_SEQS[i]
      i += 1
    return result

  def gaps(self):
    result = list(range(0,RANDOM_SEQS_LENGTH))
    i = RANDOM_SEQS_LENGTH
    for allele in self._alleles:
      i += len(allele)
      result += list(range(i,i+RANDOM_SEQS_LENGTH))
      i += RANDOM_SEQS_LENGTH
    
    return set(result)
      
  def concatenated(self):
    result = Seq("")
    for allele in self._alleles:
      result += allele
    return result
