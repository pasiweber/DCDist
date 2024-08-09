import pybind11
import numpy as np
import numba
import dctree #The homemade extension
import efficientdcdist.dctree as dcdist_pascal
from hdbscan import HDBSCAN as HDBSCAN
from HDBSCAN_nary import HDBSCAN as HDBSCANnary
from visualization import visualize, plot_tree, plot_tree_v2, plot_embedding
from kcentroids_nary import DCKCentroids as DCKCentroids_nary


import time
from benchmark import create_dataset

k = 3
dim = 2 #124 is max dimension.......
n = 20
dataset = "blobs"
points, _ = create_dataset(num_points=n, datatype=dataset, num_features=dim)

# points = np.array([   [1,2],
#                        [1,4],
#                        [2,3],
#                        [1,1],
#                        [-5,15], #5
#                        [11,13],
#                        [13,11],
#                        [10,8],
#                        [14,13],
#                        [16,17], #10
#                        [18,19],
#                        [19,18]], dtype=np.float64)
print("Testing dc_dist code with dataset", dataset, "and n =", n)


t1 = time.time()
print("1")
#Sklearn
sklearn_hdbscan = HDBSCAN(min_cluster_size=k, min_samples=k)
sklearn_hdbscan.fit(points)
sklearn_hdb_labels = sklearn_hdbscan.labels_

t2 = time.time()
print("2")
# hdbscan_nary = HDBSCANnary(min_pts=k, min_cluster_size=k, allow_single_cluster=False)
# hdbscan_nary.fit(points)
# hdb_labels = hdbscan_nary.labels_

real_k = 5
res = dctree.compute_clustering(points, k, k, real_k, ["kmedian", "kmeans", "hdbscan"])
print("fast kmedian labels:", res)
t3 = time.time()
print("3")

kmedian = DCKCentroids_nary(k=real_k, min_pts=k, loss="kmedian", noise_mode="none") #This implicitly creates the dc-tree first
kmedian.fit(points)
kmedian_labels = kmedian.labels_
print("slow kmedian labels:", kmedian_labels)
t4 = time.time()
print("4")

t5 = time.time()

# kcentroid_hierarchy = kmedian.define_cluster_hierarchy_nary(points)
# plot_tree_v2(kcentroid_hierarchy)

# plot_embedding(
#         points,
#         [hdb_labels, np.array(labels), sklearn_hdb_labels],
#         ["python slow", "c++ new", "sklearn"],
#         centers=None,
#         dot_scale=2, #6 used for small dataset, 2 used for 400 points
#         annotations=False
#     )


print("sklearn:", t2-t1)
print("old python:", t3-t2)
print("new c++", t4-t3)
#print("pargeo vector:", t5-t4)

