import read as r
import matplotlib.pyplot as plt

def cov_plot(out_dir, width=14, height=7):
  stats = r.read(out_dir + "/instance.stats.json", r.json)
  cov_sum = stats["cov_sum"]
  cov_min = stats["cov_min"]
  cov_max = stats["cov_max"]
  cov_all = stats["cov_all"]
  
  data = []
  lo = float("inf")
  hi = 0
  for k in cov_sum.keys():
    if cov_sum[k] > 0:
      if cov_min[k] < lo: lo = cov_min[k]
      if cov_max[k] > hi: hi = cov_max[k]
      data.append((k,cov_all[k]))
      
  data.sort(key = lambda a: int(a[0].split("_")[1]))
  
  plt.figure(figsize = (width, height))
  plt.boxplot([a[1] for a in data])
  plt.xticks(range(1,len(data)+1), [a[0] for a in data])
  plt.yticks(range(lo, hi+1, 20))
  plt.show()
