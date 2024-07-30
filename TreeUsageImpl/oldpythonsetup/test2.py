import pybind11
import numpy as np
import numba
import dctree #The homemade extension
import efficientdcdist.dctree as dcdist_pascal
from hdbscan import HDBSCAN as HDBSCAN
from HDBSCAN_nary import HDBSCAN as HDBSCANnary

import time
from benchmark import create_dataset



k = 3
dim = 2
n = 10
dataset = "blobs"
points1, _ = create_dataset(num_points=n, datatype=dataset, num_features=dim)
points1t = np.transpose(points1)



#dc_constructor = dcdist_pascal.DCTree(points1, k)


print("Testing dc_dist code with dataset", dataset, "and n =", n)

t1 = time.time()
#Sklearn
sklearn_hdbscan = HDBSCAN(min_cluster_size=k, min_samples=k)
sklearn_hdbscan.fit(points1)
sklearn_hdb_labels = sklearn_hdbscan.labels_

t2 = time.time()

#HDBSCAN old python
hdbscan_nary = HDBSCANnary(min_pts=k, min_cluster_size=k, allow_single_cluster=False)
hdbscan_nary.fit(points1)
hdb_labels = hdbscan_nary.labels_

t3 = time.time()

cdists3 = dctree.compute_hdbscan_labels(points1, k, "pargeo")

t4 = time.time()
cdists4 = dctree.compute_hdbscan_labels(points1, k, "arma")

t5 = time.time()




print("sklearn:", t2-t1)
print("old slow:", t3-t2)
print("pargeo orig", t4-t3)
print("pargeo arma:", t5-t4)

