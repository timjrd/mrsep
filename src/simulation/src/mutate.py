import os
import copy
import random
import subprocess

from Bio.Seq import Seq

def mutate(seq):
  r = random.randint(1,3)
  msbar = os.environ["EMBOSS_MSBAR"]
  out = subprocess.check_output([
    msbar, "-auto", "-filter",
    "-point", "4" if r == 1 else "0",
    "-block", "4" if r == 2 else "0",
    "-codon", "4" if r == 3 else "0"
  ], input=">_\n" + str(seq), encoding="ascii")
  clean_out = out.splitlines()[1:]
  result = Seq("")
  for line in clean_out: result += Seq(line)
  return result

def mutate_forward(seqs, excluded=[]):
  mut = mutate(seqs[0])
  while mut in seqs or mut in excluded:
    mut = mutate(seqs[0])
  return [mut] + seqs

# a___b
def mutate_1(seq):
  a = [seq]
  b = mutate_forward(a)
  return [b]

# a___b
#  \__c
def mutate_v(seq):
  a = [seq]
  b = mutate_forward(a)
  c = mutate_forward(a, b)
  return [b,c]

# a___b___c
def mutate_2(seq):
  a = [seq]
  b = mutate_forward(a)
  c = mutate_forward(b)
  return [c]

# a___b___d
#  \__c___e
def mutate_u(seq):
  a = [seq]
  b = mutate_forward(a)
  c = mutate_forward(a, b)
  d = mutate_forward(b, c)
  e = mutate_forward(c, d)
  return [d,e]

def mutate_0(seq):
  return [[seq]]

PATTERNS = [
  # 1 root
  ("mut_1", [mutate_1]),
  ("mut_v", [mutate_v]),
  ("mut_2", [mutate_2]),
  ("mut_u", [mutate_u]),

  # 2 roots
  ("mut_1_0", [mutate_1, mutate_0]),
  ("mut_1_1", [mutate_1, mutate_1]),
  ("mut_2_1", [mutate_2, mutate_1]),
  ("mut_v_1", [mutate_v, mutate_1]),
  ("mut_2_2", [mutate_2, mutate_2]),
  ("mut_v_2", [mutate_v, mutate_2]),
  ("mut_v_0", [mutate_v, mutate_0])
]

# seq = Seq("AAAAAAGAATTCATTATACATGATAGTTTAGTATTTGATTTGATATTAAATATAAAATTATTAAAATTCAATTTACTTGCCAATAGAAGTACTATTGGCATATTTGCTTTTATTGGTGC")
# print(seq)
# print(mutate(seq))

# max_d = 0
# for _ in range(1000):
#   seq = Seq("AAATTTCCCGGG")
#   mut = mutate(seq)
#   d = abs(len(mut)-len(seq))
#   if d > max_d: max_d = d

# print(max_d)  
