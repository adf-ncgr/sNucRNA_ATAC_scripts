#!/usr/bin/env python3
import sys
import numpy as np
import doubletdetection
import tarfile
import matplotlib.pyplot as plt

matrix_path = sys.argv[1]
raw_counts = doubletdetection.load_mtx(matrix_path)
zero_genes = np.sum(raw_counts, axis=0) == 0
raw_counts = raw_counts[:, ~zero_genes]
clf = doubletdetection.BoostClassifier(n_iters=50)
doublets = clf.fit(raw_counts).predict(p_thresh=1e-7, voter_thresh=0.8)

db = doublets.nonzero()
f = open(sys.argv[2], "w")
for d in db[0] :
    f.write(str(d+1) + "\n")
f.close()
