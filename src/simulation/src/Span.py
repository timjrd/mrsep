
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
  
  def from_record(self, strain, aligned, record):
    self.clear()
    
    interleaved = strain.interleaved()
    gaps        = strain.gaps()
    
    start_il = None
    if record.is_reversed:
      start_il = len(interleaved) - len(record.read) - record.offset
    else:
      start_il = record.offset
    
    pos_al = -1
    for pos_il in range(len(interleaved)):
      if pos_il in gaps: pass
      else:
        pos_al += 1
        if pos_al >= len(aligned): break
        # print(str(pos_al) + "/" + str(len(aligned)))
        while aligned[pos_al] == "-":
          pos_al += 1
        if pos_il >= start_il and pos_il < start_il + len(record.read):
          self.append(self.Pos(pos_al, record.read[pos_il - start_il], aligned[pos_al]))
      
    return self
