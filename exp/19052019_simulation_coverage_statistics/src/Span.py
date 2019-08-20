
class Span(list):
  class Pos:    
    def __init__(self, pos=None, read_base=None, ref_base=None):
      self.pos       = pos
      self.read_base = read_base
      self.ref_base  = ref_base

    def to_text(self):
      return str(self.read_base) + " " + str(self.pos + 1)

    def from_json(self, data):
      vars(self).update(data)
      return self
  
  def to_text(self):
    return " ".join(map(lambda x: x.to_text(), self))

  def to_json(self):
    return list(map(vars,self))

  def from_json(self, data):
    self.clear()
    self.extend(map(lambda x: self.Pos().from_json(x), data))
    return self
  
  def from_record(self, params, record):
    self.clear()
    
    offset = record.offset - params.ends_length
    edge   = offset < 0 and abs(offset) < len(record.read)
    
    read = None
    if not edge:
      read = record.read
    elif record.is_reversed:
      read = record.read[:offset]
    else:
      read = record.read[abs(offset):]
    
    offset = 0 if edge else offset
    
    if offset >= 0:
      ref = list(filter(
        lambda x: x[1] != "-",
        enumerate(params.aligned_alleles()[record.allele_id].seq)
      ))
      
      offset   = len(ref)-len(read)-offset if record.is_reversed else offset
      overflow = len(ref)-len(read)-offset
      
      if offset < 0:
        read = read[abs(offset):]
      elif overflow < 0:
        read = read[:overflow]
        
      offset = max(0,offset)
      
      for i,x in enumerate(read):
        r = ref[offset+i]
        self.append(self.Pos(r[0], x, r[1]))
    
    return self
