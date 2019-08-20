import sys
import random
import copy
import json

from Bio.Seq import Seq

def _random_seq(length):
  return Seq("".join([
    random.choice("ATCG") for _ in range(length)
  ]))

class Sample(list):
  def new(self, params):
    self.clear()
    
    random.seed(params.seed)
    my_alleles = copy.deepcopy(params.alleles())
    random.shuffle(my_alleles)
    chosen = my_alleles[0 : params.nb_alleles]
    
    for a in chosen:
      n = random.randint(params.min_coby_nb, params.max_copy_nb)
      for i in range(n):
        b = copy.deepcopy(a)
        l = _random_seq(params.ends_length)
        r = _random_seq(params.ends_length)
        b.seq = l + b.seq + r
        self.append(b)
    
    random.shuffle(self)
    return self
  
  def from_fasta(self, data):
    self.clear()
    self.extend(data)
    return self
