from multiprocessing import Pool

import read as r
from Params import Params
from Simulation import Simulation

seeds = [21249066, 86405941]

def job(gene):
  for take_ref in [True,False]:
    for seed in seeds:
      r = "ref" if take_ref else "sim"
      
      p = Params()
      p.alleles_file         = "data/" + gene + ".fasta"
      p.aligned_alleles_file = "data/" + gene + "_mafft.fasta"
      p.nb_alleles  = 10
      p.min_coby_nb = 1
      p.max_copy_nb = 3
      p.ends_length = 200
      p.seed        = seed
      p.art_flags   = ["--seqSys","HS20","--len","100","--fcov","100"]
      p.out_dir     = f"results/{gene}/" + p.hash() + f"/{r}"
      p.take_ref    = take_ref
      
      Simulation().new(p)

genes = r.read("data/genes.txt", r.text).splitlines()
genes = map(lambda x: x.strip(), genes)

with Pool() as pool:
  pool.map(job, genes)
