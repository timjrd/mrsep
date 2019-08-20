import os
import subprocess

import samples as s
import read  as r
import write as w
from Params import Params
from Sample import Sample
from Instance import Instance

class Simulation:
  def new(self, root_dir, strains, aligned, sample):
    out_dir = root_dir + "/" + sample[0]
    
    os.makedirs(out_dir, exist_ok=True)
    
    # sample_json = out_dir + "/sample.json"
    # w.write(sample_json, w.json(s.to_json(sample)))

    sample_fasta = out_dir + "/sample.fasta"
    w.write(sample_fasta, w.fasta(s.to_fasta(sample)))

    art_prefix = out_dir + "/art"
    art = os.environ['ART_ILLUMINA']
    subprocess.run([art,
      "--in"    , sample_fasta,
      "--out"   , art_prefix,
                    
      "--seqSys", "HS20",
      "--len"   , "100",
      "--fcov"  , "100"                    
    ], stdout=subprocess.DEVNULL)
    take_ref   = False
    art_output = r.read(art_prefix + ".aln", r.aln(take_ref))

    instance = Instance().new(strains, aligned, art_output)
    # w.write(out_dir + "/instance.json", w.json(instance.to_json()))
    w.write(out_dir + "/instance.txt", w.text(instance.to_text()))
    w.write(out_dir + "/instance.stats.json", w.json(instance.stats()))
  
  # def load(self, out_dir):
  #   self.params = Params().from_json(
  #     r.read(out_dir + "/params.json", r.json)
  #   )
  #   self.sample = Sample().from_fasta(
  #     r.read(out_dir + "/sample.fasta", r.fasta_list)
  #   )
  #   self.art_output = r.read(out_dir + "/art.aln", r.aln(self.params.take_ref))
  #   self.instance = Instance().from_json(
  #     r.read(out_dir + "/instance.json", r.json)
  #   )
  #   return self
