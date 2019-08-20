import os.path

import read as r
from Strain import Strain

def base_strains(directory):
  result      = []
  
  genes       = []
  genes_names = []
  
  strains_db = r.read(directory + "/ST.txt", r.csv)
  
  for gene in strains_db[0]:
    gene_file = directory + "/" + gene + ".fasta"
    if os.path.isfile(gene_file):
      genes.append( r.read(gene_file, r.fasta_dict) )
      genes_names.append(gene)

  for strain in strains_db[1:]:
    alleles = []
    id = strain[0]
    gene = 0
    
    for allele in strain[1:]:
      try:
        name = genes_names[gene] + "_" + allele
        alleles.append(genes[gene][name].seq)
      except KeyError  : pass
      except IndexError: pass
      gene += 1
    
    result.append(Strain(id,alleles))

  return result
