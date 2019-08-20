import json as j
import csv as c

from Bio.Seq import Seq
from Bio import SeqIO

def read(filename, reader):
  r = None
  with open(filename) as f:
    r = reader(f)
  return r

def json(f):
  return j.load(f)

def fasta_list(f):
  return list(SeqIO.parse(f,"fasta"))

def fasta_dict(f):
  return SeqIO.to_dict(SeqIO.parse(f,"fasta"))

def text(f):
  return f.read()

def csv(f):
  result = []
  r = c.reader(f, delimiter="\t")
  for row in r:
    row_ = []
    for col in row:
      row_.append(col)
    result.append(row_)
    
  return result

class AlnRecord:
  def __init__(self):
    self.strain_id   = None
    self.offset      = None
    self.is_reversed = None
    self.read        = None
  
def aln(take_ref):
  def reader(f):
    result = []
    record = AlnRecord()
    skip   = None
    for line in f:
      if line[0] == ">":
        words = line[1:].split()
        record.strain_id = words[0]
        record.offset = int(words[2])
        record.is_reversed = words[3] == "-"
        skip = not take_ref
      elif skip == True:
        skip = False
      elif skip == False:
        read = Seq(line.strip())
        record.read = read.reverse_complement() if record.is_reversed else read
        result.append(record)
        record = AlnRecord()
        skip   = None
        
    return result

  return reader
