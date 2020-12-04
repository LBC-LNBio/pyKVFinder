import os
import sys
import time
import numpy

# sys.argv[1]: kv1000 directory
# sys.argv[2]: output directory
# sys.argv[3]: number of threads

pdbs = [pdb for pdb in sorted(os.listdir(f"{sys.argv[1]}")) if pdb.endswith('.pdb')]
# print(pdbs)

fn = os.path.join(sys.argv[2], f"{int(sys.argv[3]):02d}.csv")
with open(fn, 'w') as f:
    f.write("pdb, avg_time, std_time\n")
    for pdb in pdbs:
        print(pdb)
        t = numpy.zeros(5)
        for i in range(0, 5):
            start = time.time()
            os.system(f"pyKVFinder {os.path.join(sys.argv[1], pdb)} --nthreads {sys.argv[3]}")
            t[i] = time.time() - start
        f.write(f"{pdb.replace('.pdb', '')}, {t.mean():.4f}, {t.std():.4f}\n")
