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
dim = 200 #124 is max dimension.......
n = 1000
dataset = "blobs"
points1, _ = create_dataset(num_points=n, datatype=dataset, num_features=dim)
points1t = np.transpose(points1)

points2 = np.array([   [1.1, 2.5],
                       [1.4, 4.4],
                       [2.2, 3.4],
                       [1.2, 1.3],
                       [-5.4,15.1], #5
                       [11.6,13.5],
                       [13.5,11.6],
                       [10.2,8.5],
                       [14.1,13.5],
                       [16.4,17.4], #10
                       [18.2,19.1],
                       [19.4,18.3]], dtype=np.float64) #If points are given as ints, we need to explicitly convert them to float.
points2t = np.transpose(points2)

#dc_constructor = dcdist_pascal.DCTree(points1, k)


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

cdists3 = dctree.compute_hdbscan_labels(points1, k, "arma")


t4 = time.time()
print("4")
#cdists4 = dctree.compute_hdbscan_labels(points1t, k, "vector")

t5 = time.time()




print("sklearn:", t2-t1)
print("pargeo:", t3-t2)
print("pargeo arma", t4-t3)
#print("pargeo vector:", t5-t4)

