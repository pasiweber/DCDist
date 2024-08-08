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
dim = 20 #124 is max dimension.......
n = 10000
dataset = "blobs"
points1, _ = create_dataset(num_points=n, datatype=dataset, num_features=dim)


points2 = np.array([[1,2],
                       [1,4],
                       [2,3],
                       [1,1],
                       [-5,15], #5
                       [11,13],
                       [13,11],
                       [10,8],
                       [14,13],
                       [16,17], #10
                       [18,19],
                       [19,18]], dtype=np.float64)



print("Testing dc_dist code with dataset", dataset, "and n =", n)

t1 = time.time()
print("1")
#Sklearn
sklearn_hdbscan = HDBSCAN(min_cluster_size=k, min_samples=k)
sklearn_hdbscan.fit(points1)
sklearn_hdb_labels = sklearn_hdbscan.labels_

t2 = time.time()
print("2")
#cdists2 = dctree.compute_hdbscan_labels(points1, k, "pargeo")

t3 = time.time()
print("3")

labels = dctree.compute_hdbscan_labels(points1, k, k)


t4 = time.time()
print("4")
#cdists4 = dctree.compute_hdbscan_labels(points1t, k, "vector")

t5 = time.time()




print("sklearn:", t2-t1)
print("pargeo:", t3-t2)
print("pargeo arma", t4-t3)
#print("pargeo vector:", t5-t4)

