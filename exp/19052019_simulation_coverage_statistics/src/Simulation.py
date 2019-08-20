import os
import subprocess

import read  as r
import write as w
from Params import Params
from Sample import Sample
from Instance import Instance

class Simulation:
  def new(self, params):
    os.makedirs(params.out_dir, exist_ok=True)
    
    self.params = params
    w.write(params.out_dir + "/params.json", w.pretty_json(params.to_json()))
    
    self.sample = Sample().new(params)
    sample_file = params.out_dir + "/sample.fasta"
    w.write(sample_file, w.fasta(self.sample))
    
    art_prefix = params.out_dir + "/art"
    art = os.environ['ART_ILLUMINA']
    subprocess.run([art,
      "--in"     , sample_file,
      "--out"    , art_prefix,
      "--rndSeed", str(params.seed)
    ] + params.art_flags, stdout=subprocess.DEVNULL)
    self.art_output = r.read(art_prefix + ".aln", r.aln(params.take_ref))
    
    self.instance = Instance().new(params, self.art_output)
    w.write(params.out_dir + "/instance.json", w.json(self.instance.to_json()))
    w.write(params.out_dir + "/instance.txt", w.text(self.instance.to_text()))
    w.write(params.out_dir + "/instance.stats.json", w.json(self.instance.stats()))

    return self
  
  def load(self, out_dir):
    self.params = Params().from_json(
      r.read(out_dir + "/params.json", r.json)
    )
    self.sample = Sample().from_fasta(
      r.read(out_dir + "/sample.fasta", r.fasta_list)
    )
    self.art_output = r.read(out_dir + "/art.aln", r.aln(self.params.take_ref))
    self.instance = Instance().from_json(
      r.read(out_dir + "/instance.json", r.json)
    )
    return self
