import os
import sys

import read as r
import write as w
import base_strains as bs
import samples as s
import alignment as al

import Strain
from Simulation import Simulation

from Bio.Seq import Seq

root_dir     = sys.argv[1]
strains_file = root_dir + "/strains.json"
aligned_file = root_dir + "/aligned.json"
samples_file = root_dir + "/samples.json"

def ko(path):
  return not os.path.isfile(path)

if ko(strains_file) or ko(aligned_file) or ko(samples_file):
  base_strains = bs.base_strains("../../data/Borrelia/MLST_19032019")
  # base_strains = base_strains[:22]
  
  samples = s.samples(base_strains)

  more_strains = {}
  for sample in samples:
    for strain in sample[1]:
      more_strains[strain.id] = strain
  more_strains = list(more_strains.values())

  strains = base_strains + more_strains
  aligned = al.align(strains)
  
  os.makedirs(root_dir, exist_ok=True)
  w.write(strains_file, w.json(list(map(lambda s: s.to_json(), strains))))
  w.write(aligned_file, w.json({k : str(v) for k,v in aligned.items()}))
  w.write(samples_file, w.json(list(map(s.to_json, samples))))

strains_json = r.read(strains_file, r.json)
aligned_json = r.read(aligned_file, r.json)
samples_json = r.read(samples_file, r.json)

strains = list(map(Strain.from_json, strains_json))
aligned = {k : Seq(v) for k,v in aligned_json.items()}
samples = list(map(s.from_json, samples_json))

for sample in samples:
  Simulation().new(root_dir, strains, aligned, sample)
