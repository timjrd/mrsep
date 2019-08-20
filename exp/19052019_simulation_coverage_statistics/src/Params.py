import json
import hashlib

import read as r

class Params:
  def __init__(self):
    self.alleles_file         = None
    self.aligned_alleles_file = None
    self.out_dir              = None
    
    self.nb_alleles  = None
    self.min_coby_nb = None
    self.max_copy_nb = None
    self.ends_length = None
    self.take_ref    = None
    self.seed        = None
    self.art_flags   = None
    
    self._alleles         = None
    self._aligned_alleles = None
    
  def alleles(self):
    if self._alleles == None and self.alleles_file != None:
      self._alleles = r.read(self.alleles_file, r.fasta_list)
        
    return self._alleles

  def aligned_alleles(self):
    if self._aligned_alleles == None and self.aligned_alleles_file != None:
      self._aligned_alleles = r.read(self.aligned_alleles_file, r.fasta_dict)
    
    return self._aligned_alleles
  
  def to_json(self):
    return { k:v for k,v in vars(self).items() if not k.startswith("_") }
  
  def from_json(self, data):
    vars(self).update(data)
    return self
    
  def hash(self):
    s = json.dumps(self.to_json(), sort_keys=True)
    return hashlib.sha256(s.encode("utf-8")).hexdigest()[:10]
