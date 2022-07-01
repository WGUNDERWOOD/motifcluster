from matplotlib import rc
from matplotlib import pylab
import seaborn
import numpy as np
import pandas as pd

# to fix fonts:
# may need to use either
# 'CMU Serif' or 'Computer Modern'
# depending on operating system.
###########################################
rc('text', usetex=True)
seaborn.set()
rc('font', **{'family': 'serif', 'serif': ['CMU Serif']})
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
###########################################

ks = [10, 100]
langs = ["python", "r"]
graph_types = ["erdos_renyi", "barabasi_albert"]

for graph_type in graph_types:
  for lang in langs:
    for k in ks:

      # plot graph
      pylab.figure(figsize=(10,4))

      for method in ["dense", "sparse"]:

        # read data
        results_filename = "results/" + lang + "_k" + str(k) + "_" + method + "_" + graph_type + ".csv"
        results = pd.read_csv(results_filename)

        for motif in ["M1", "M8", "M11"]:

          # plot line
          line_results = results[results["motif"] == motif]
          line_results = line_results[line_results["method"] == method]
          plot_label = "$\mathcal{M}_{" + motif[1:] + "}$"

          xs = []
          ys = []
          es = []
          for n in sorted(set(line_results["n"])):
            point_results = line_results[line_results["n"] == n]
            xs.append(n)
            ys.append(np.mean(point_results["time"]))
            es.append(np.std(point_results["time"]))

          pylab.errorbar(xs, ys, es, label=plot_label)

      pylab.legend(ncol=2, title='Dense  \\hspace{1.2cm}        Sparse     \\hspace{0.4cm}  ')
      pylab.xlabel('$n$', fontsize=18)
      pylab.ylabel('Time (s)', fontsize=18)
      pylab.xscale('log')
      pylab.yscale('log')

      if lang == "python":
        if graph_type == "erdos_renyi":
          if k == 10:
            pylab.title('Python timing results for the sparser regime $\mathrm{ER}(n,' + str(k) + '/n)$',fontsize=18)
          elif k == 100:
            pylab.title('Python timing results for the less sparse regime $\mathrm{ER}(n,' + str(k) + '/n)$',fontsize=18)
        elif graph_type == "barabasi_albert":
          if k == 10:
            pylab.title('Python timing results for the sparser regime $\mathrm{BA}(n,' + str(k) + ')$',fontsize=18)
          elif k == 100:
            pylab.title('Python timing results for the less sparse regime $\mathrm{BA}(n,' + str(k) + ')$',fontsize=18)

      if lang == "r":
        if graph_type == "erdos_renyi":
          if k == 10:
            pylab.title('R timing results for the sparser regime $\mathrm{ER}(n,' + str(k) + '/n)$',fontsize=18)
          elif k == 100:
            pylab.title('R timing results for the less sparse regime $\mathrm{ER}(n,' + str(k) + '/n)$',fontsize=18)
        elif graph_type == "barabasi_albert":
          if k == 10:
            pylab.title('R timing results for the sparser regime $\mathrm{BA}(n,' + str(k) + ')$',fontsize=18)
          elif k == 100:
            pylab.title('R timing results for the less sparse regime $\mathrm{BA}(n,' + str(k) + ')$',fontsize=18)

      pylab.xticks(fontsize=16)
      pylab.yticks(fontsize=16)
      plot_filename = "plots/" + lang + "_k" + str(k) + "_" + graph_type + ".pdf"
      pylab.savefig(plot_filename, bbox_inches="tight")
      pylab.close('all')
