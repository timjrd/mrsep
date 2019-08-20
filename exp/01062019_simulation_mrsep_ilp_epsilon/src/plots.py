import json
import os
import subprocess

import matplotlib.pyplot as plt


WIDTH=20
HEIGHT=13


def plot_all(root, width=WIDTH, height=HEIGHT):
  out = subprocess.check_output(
    "find {} -name solver.json".format(root),
    stderr=subprocess.DEVNULL,
    shell=True,
    universal_newlines=True
  )
  
  paths = [x[len(root)+1:].split("/") for x in out.splitlines()]
  for x in paths:
    x.pop()
    x.pop(2)
  
  results = []
  for x in paths:
    if x not in results:
      results.append(x)

  results.sort()

  for x in results:
    prefix = "/".join(x[:2])
    suffix = "/".join(x[-1:])
    ref_path = "{}/{}/ref/{}".format(root,prefix,suffix)
    sim_path = "{}/{}/sim/{}".format(root,prefix,suffix)
    plot(root, ref_path, sim_path, width, height)


def frange(start, end, step):
  x = start
  xs = []
  while x < end:
    xs.append(x)
    x += step    
  return xs


def plot(root, ref_path, sim_path, width=WIDTH, height=HEIGHT):
  plt.rcParams.update({
    "font.size": 22
  })
  fig, ax = plt.subplots(figsize = (width,height))
  
  p1 = plot1(ax, ref_path, "reference", 10, width, height)
  p2 = plot1(ax, sim_path, "simulated", 25, width, height)
  
  params = None
  stats = None
  ref_params = None
  sim_params = None
  ref_solver = None  
  sim_solver = None
  with open(root     + "/params.json")         as f: params     = json.load(f)
  with open(sim_path + "/instance.stats.json") as f: stats      = json.load(f)
  with open(ref_path + "/params.json")         as f: ref_params = json.load(f)
  with open(sim_path + "/params.json")         as f: sim_params = json.load(f)
  with open(ref_path + "/solver.json")         as f: ref_solver = json.load(f)
  with open(sim_path + "/solver.json")         as f: sim_solver = json.load(f)
  
  print()
  print("Paths              {}".format(ref_path.replace("ref", "*")))
  print("Datasets           {} {}".format(ref_params["dataset"], sim_params["dataset"]))
  print("Alleles database   {}".format(ref_params["aligned_nb_alleles"]))
  print("Alleles in reads   {}".format(ref_params["reads_nb_alleles"]))
  print("Reads              {}".format(stats["nb_reads"]))
  print("Bases              {}".format(stats["nb_bases"]))
  print("Faulty reads       {:.2f} %".format(100 * float(stats["faulty_reads"]) / float(stats["nb_reads"])))
  print("Faulty bases       {:.2f} %".format(100 * float(stats["faulty_bases"]) / float(stats["nb_bases"])))
  print("Timeout            {} minutes".format(params["timeout"]))
  
  print()
  print("Solver stats (sim):")
  stats = {}
  for e in sim_solver:
    if "stats" in e:
      stats = e["stats"]
      break
  for k,v in stats.items():
    print("  {: <6} {}".format(v,k))
  
  print()
  print("Solver stats (ref):")
  stats = {}
  for e in ref_solver:
    if "stats" in e:
      stats = e["stats"]
      break
  for k,v in stats.items():
    print("  {: <6} {}".format(v,k))
  
  if p1 == None and p2 == None:
    plt.close(fig)
    print()
    print("NO RESULTS")
  else:
    n  = (float("inf"), -1, float("inf"), -1)
    p1 = n if p1 == None else p1
    p2 = n if p2 == None else p2
    
    x_lo = min(p1[0], p2[0])
    x_hi = max(p1[1], p2[1])
    y_lo = min(p1[2], p2[2])
    y_hi = max(p1[3], p2[3])
    
    m = max(1, int((y_hi - y_lo) / 10))
    ax.set_xticks(frange(x_lo, x_hi+0.06, 0.05))
    ax.set_yticks(range(y_lo, y_hi+m, m))
    
    ax.set_xlabel("epsilon"    , labelpad=30)
    ax.set_ylabel("assignments", labelpad=30)
    ax.legend(loc="lower right")
    
    plt.show()
    plt.close(fig)
    
  print()
  print()


def plot1(ax, path, label, linewidth, width=WIDTH, height=HEIGHT):
  data = None
  with open(path + "/solver.json") as f: data = json.load(f)
  
  xs = [    d["epsilon" ]  for d in data if "solution" in d]
  ys = [len(d["solution"]) for d in data if "solution" in d]
  if xs and ys:
    ax.scatter(xs, ys, label=label, linewidth=linewidth)
    x_lo = min(xs)
    x_hi = max(xs)
    y_lo = min(ys)
    y_hi = max(ys)
    return (x_lo, x_hi, y_lo, y_hi)
  else:
    return None
