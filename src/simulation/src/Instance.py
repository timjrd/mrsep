from collections import defaultdict
import random

from Bio.Seq import Seq

from Span   import Span
from Strain import Strain

class Instance:
  def new(self, strains, aligned, records):
    self.strains = { s.id : s for s in strains }
    self.aligned = aligned
    self.reads   = defaultdict(lambda: [])
    
    for record in records:
      span = Span().from_record(
        self.strains[record.strain_id],
        self.aligned[record.strain_id],
        record
      )
      if span:
        self.reads[record.strain_id].append(span)
    
    return self
  
  def filter(self, strains_nb_strains, reads_nb_strains):
    result = Instance()
    
    strains_ids = list(self.reads.keys())
    strains_ids.sort()
    random.shuffle(strains_ids)
    
    more_ids = [ k for k in self.strains.keys() if k not in strains_ids ]
    random.shuffle(more_ids)
    strains_ids.extend(more_ids)
    
    strains_chosen = strains_ids[0 : strains_nb_strains]
    reads_chosen   = strains_ids[0 : reads_nb_strains  ]
    
    result.strains = { k:v for k,v in self.strains.items() if k in strains_chosen }
    result.aligned = { k:v for k,v in self.aligned.items() if k in strains_chosen }
    result.reads   = { k:v for k,v in self.reads.items()   if k in reads_chosen   }
    
    return result
  
  # def to_json(self):
  #   return {
  #     "strains": { k : v.to_json() for k,v in self.strains.items() },
  #     "reads":   { k : list(map(lambda x: x.to_json(), v)) for k,v in self.reads.items() }
  #   }

  # Ignoring restoration for now
  # def from_json(self, data):
  #   self.aligned_alleles = { k : Strain(v) for k,v in data["aligned_alleles"].items() }
  #   self.reads = { k : [Span().from_json(s) for s in v] for k,v in data["reads"].items() }
  #   return self
      
  def to_text(self):
    all_reads = []
    keys = list(self.reads.keys())
    keys.sort()
    for k in keys:
      all_reads.extend(self.reads[k])
    
    all_aligned = []
    keys = list(self.aligned.keys())
    keys.sort()
    for k in keys:
      all_aligned.append(self.aligned[k])

    return "\n".join(
      list(map(str, all_aligned)) + [""] +
      list(map( lambda ix: str(ix[0] + 1) + " 1 " + ix[1].to_text(),
                enumerate(all_reads) ))
    )
  
  def stats(self):
    result = defaultdict(lambda: {})
    
    cov = defaultdict(lambda: 0)
    faulty_bases = 0
    faulty_reads = 0
    max_faulty_bases_per_read = 0
    for strain_id in self.reads.keys():
      for read in self.reads[strain_id]:
        faulty_bases_in_read = 0
        for x in read:
          if x.read_base == x.ref_base:
            cov[(strain_id, x.pos)] += 1
          else:
            faulty_bases_in_read += 1
        if faulty_bases_in_read > 0:
          faulty_bases += faulty_bases_in_read
          faulty_reads += 1
          if faulty_bases_in_read > max_faulty_bases_per_read:
            max_faulty_bases_per_read = faulty_bases_in_read
    
    result["faulty_bases"             ] = faulty_bases
    result["nb_bases"                 ] = sum(map(lambda x: sum(map(len,x)), self.reads.values()))
    result["faulty_reads"             ] = faulty_reads
    result["nb_reads"                 ] = sum(map(len,self.reads.values()))
    result["max_faulty_bases_per_read"] = max_faulty_bases_per_read
    
    for strain_id, seq in self.aligned.items():
      lo = float("inf")
      hi = 0
      s  = 0
      n  = 0
      cs = []
      for pos, base in enumerate(seq):
        if base != "-":
          c = cov[(strain_id, pos)]
          cs.append(c)
          if c < lo: lo = c
          if c > hi: hi = c
          s += c
          n += 1
          
      result["cov_min"][strain_id] = lo
      result["cov_max"][strain_id] = hi
      result["cov_sum"][strain_id] = s
      result["cov_nb" ][strain_id] = n
      result["cov_all"][strain_id] = cs
    
    return result
