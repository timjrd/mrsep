from collections import defaultdict
import random

from Bio.Seq import Seq

from Span import Span

class Instance:
  def new(self, params, art_output):
    self.aligned_alleles = { k:v.seq for k,v in params.aligned_alleles().items() }
    self.reads = defaultdict(lambda: [])
    
    for record in art_output:
      span = Span().from_record(params, record)
      if span: self.reads[record.allele_id].append(span)
    
    return self
  
  def filter(self, params, aligned_nb_alleles, reads_nb_alleles):
    random.seed(params.seed)
    
    result = Instance()
    
    allele_ids = list(self.reads.keys())
    allele_ids.sort()
    random.shuffle(allele_ids)
    
    more_ids = [ k for k in self.aligned_alleles.keys() if k not in allele_ids ]
    random.shuffle(more_ids)
    allele_ids.extend(more_ids)
    
    aligned_chosen = allele_ids[0 : aligned_nb_alleles]
    reads_chosen   = allele_ids[0 : reads_nb_alleles  ]
    
    result.aligned_alleles = { k:v for k,v in self.aligned_alleles.items() if k in aligned_chosen }
    result.reads           = { k:v for k,v in self.reads.items()           if k in reads_chosen   }
    
    return result
  
  def to_json(self):
    return {
      "aligned_alleles": { k : str(v) for k,v in self.aligned_alleles.items() },
      "reads": { k : list(map(lambda x: x.to_json(), v)) for k,v in self.reads.items() }
    }

  def from_json(self, data):
    self.aligned_alleles = { k : Seq(v) for k,v in data["aligned_alleles"].items() }
    self.reads = { k : [Span().from_json(s) for s in v] for k,v in data["reads"].items() }
    return self
      
  def to_text(self):
    all_reads = []
    keys = list(self.reads.keys())
    keys.sort()
    for k in keys:
      all_reads.extend(self.reads[k])
    
    aligned_alleles = []
    keys = list(self.aligned_alleles.keys())
    keys.sort()
    for k in keys:
      aligned_alleles.append(self.aligned_alleles[k])
    
    return "\n".join(
      list(map(str, aligned_alleles)) + [""] +
      list(map( lambda ix: str(ix[0] + 1) + " 1 " + ix[1].to_text(),
                enumerate(all_reads) ))
    )
  
  def stats(self):
    result = defaultdict(lambda: {})
    
    cov = defaultdict(lambda: 0)
    faulty_bases = 0
    faulty_reads = 0
    max_faulty_bases_per_read = 0
    for allele_id in self.reads.keys():
      for read in self.reads[allele_id]:
        faulty_bases_in_read = 0
        for x in read:
          if x.read_base == x.ref_base:
            cov[(allele_id, x.pos)] += 1
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
    
    for allele_id, seq in self.aligned_alleles.items():
      lo = float("inf")
      hi = 0
      s  = 0
      n  = 0
      cs = []
      for pos, base in enumerate(seq):
        if base != "-":
          c = cov[(allele_id, pos)]
          cs.append(c)
          if c < lo: lo = c
          if c > hi: hi = c
          s += c
          n += 1
          
      result["cov_min"][allele_id] = lo
      result["cov_max"][allele_id] = hi
      result["cov_sum"][allele_id] = s
      result["cov_nb" ][allele_id] = n
      result["cov_all"][allele_id] = cs
    
    return result
