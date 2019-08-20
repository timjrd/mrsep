import matplotlib.pyplot as plt
from collections import defaultdict

def plots(filename, width = 20, height = 10):
  raw   = None
  stats = defaultdict(lambda: defaultdict(lambda: {}))
  with open(filename) as f: raw = f.read()
  for line in raw.splitlines():
    ks, v = line.split()
    value = int(v)
    keys  = ks.split(".")
    k1    = keys[0]
    k2    = keys[1]
    if len(keys) == 2:
      stats[k1][k2] = value
    elif len(keys) == 3:
      if value > 0:
        k3 = int(keys[2])
        stats[k1][k2][k3] = value
    else:
      raise "malformed stats file"
  
  # p_xs = []
  # p_ys = []
  # for k,v in sorted(stats["unfiltered"]["min_perfect"].items(), key=lambda x: x[1]):
  #   p_xs.append(str(k))
  #   p_ys.append(v)
  
  # u_xs = []
  # u_ys = []
  # for k,v in sorted(stats["unfiltered"]["min_unambiguous"].items(), key=lambda x: x[1]):
  #   u_xs.append(str(k))
  #   u_ys.append(v)    
  
  # fig = plt.figure(figsize = (width,height))
  # plt.xlabel("strain type id"  , labelpad=15)
  # plt.ylabel("minimum coverage", labelpad=15)

  # plt.title("before filtering", pad=15)
  # plt.scatter(p_xs, p_ys, label="perfect matches", linewidth=9, c="black")
  # plt.scatter(u_xs, u_ys, label="perfect & unambiguous", linewidth=2, c="orange")
  # plt.xticks(p_xs)
  # ys = p_ys if p_ys else u_ys
  # plt.yticks(range(min(ys), max(ys)+25, 25))
  # plt.legend(loc="upper left")
  
  # plt.show()
  # plt.close(fig)

  p_xs = []
  p_ys = []
  for k,v in sorted(stats["filtered"]["min_perfect"].items(), key=lambda x: x[1]):
    p_xs.append(str(k))
    p_ys.append(v)
  
  u_xs = []
  u_ys = []
  for k,v in sorted(stats["filtered"]["min_unambiguous"].items(), key=lambda x: x[1]):
    u_xs.append(str(k))
    u_ys.append(v)
  
  fig = plt.figure(figsize = (width,height))
  plt.xlabel("strain type id"  , labelpad=15)
  plt.ylabel("minimum coverage", labelpad=15)

  plt.title("after filtering", pad=15)
  plt.scatter(p_xs, p_ys, label="perfect matches", linewidth=9, c="black")
  plt.scatter(u_xs, u_ys, label="perfect & unambiguous", linewidth=2, c="orange")
  plt.xticks(p_xs)
  ys = p_ys if p_ys else u_ys
  plt.yticks(range(min(ys), max(ys)+25, 25))
  plt.legend(loc="upper left")
  
  plt.show()
  plt.close(fig)
