import os
import json as j

from Bio import SeqIO

def write(filename, writer):
  r = None
  os.makedirs(os.path.dirname(filename), exist_ok=True)
  with open(filename, "w") as f:
    r = writer(f)
  return r

def json(data):
  return (lambda f: j.dump(data, f, sort_keys=True, separators=(",",":")))

def pretty_json(data):
  return (lambda f: j.dump(data, f, sort_keys=True, indent=2))

def fasta(data):
  return (lambda f: SeqIO.write(data, f, "fasta"))

def text(data):
  return (lambda f: f.write(data))
