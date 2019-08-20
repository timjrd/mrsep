from collections import defaultdict
import time
import os
import subprocess
from multiprocessing import Pool

import write as w
from Simulation import Simulation


def frange(start, end, step):
  x = start
  xs = []
  while x <= end:
    xs.append(x)
    x += step    
  return xs


####################
# BEGIN PARAMETERS #
####################
datasets           = ["clpA/3721cdda46/sim", "clpA/3721cdda46/ref"]
epsilons           = [0.1, 0.2, 0.3]
timeout            = 180
aligned_nb_alleles = [10, 50, 100]
reads_nb_alleles   = [1, 2, 3]
####################
# END   PARAMETERS #
####################


args = []
for dataset in datasets:
  sim = Simulation().load("data/" + dataset)
  for an in aligned_nb_alleles:
    for rn in reads_nb_alleles:
      rn = min(rn, len(sim.instance.reads))
      an = max(min(an, len(sim.instance.aligned_alleles)), rn)
      
      out_dir = "results/{}/a{:03}-r{:03}".format(dataset, an, rn)
      try:
        os.makedirs(out_dir)
      except OSError:
        out_dir = None
          
      if out_dir != None:
        filtered = sim.instance.filter(sim.params, an, rn)
        w.write(out_dir + "/instance.json", w.json(filtered.to_json()))
        w.write(out_dir + "/instance.stats.json", w.json(filtered.stats()))
        w.write(out_dir + "/params.json", w.pretty_json({
          "dataset"           : dataset,
          "aligned_nb_alleles": an,
          "reads_nb_alleles"  : rn
        }))
        input_file = out_dir + "/instance.txt"
        w.write(input_file, w.text(filtered.to_text()))
        for e in epsilons:
          args.append((e, input_file, out_dir))


w.write("results/params.json", w.pretty_json({
  "timeout": timeout,
  "nb_jobs": len(args)
}))


def job(args):
  start_time = time.time()
  
  epsilon, input_file, out_dir = args
  cmd = [os.environ["SOLVER"], str(epsilon), input_file]
  proc = subprocess.Popen(cmd,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    universal_newlines=True
  )
  try:
    out, err = proc.communicate(timeout=60*timeout)
    
    solution = [list(map(int,line.split())) for line in out.splitlines()[1:]]
    stats = { line.split()[0] : int(line.split()[1]) for line in err.splitlines() }
    
    solution.sort()
    
    return (out_dir, {
      "epsilon" : epsilon,
      "solution": solution,
      "stats"   : stats,
      "duration": time.time() - start_time
    })
  except subprocess.TimeoutExpired:
    proc.kill()
    return (out_dir, {
      "epsilon" : epsilon,
      "duration": time.time() - start_time
    })


def progress(i):
  print("\r {:03} % ".format(int(100.0 * float(i) / float(len(args)))), end="", flush=True)


with Pool() as pool:
  result = defaultdict(lambda: [])
  
  i = 0
  progress(i)
  for x in pool.imap_unordered(job, args):
    result[x[0]].append(x[1])    
    i += 1
    progress(i)
  
  for out_dir, r in result.items():
    w.write(out_dir + "/solver.json", w.json(r))

print()
