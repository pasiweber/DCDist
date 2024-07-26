import pybind11
import numpy as np
import numba
import dctree #The homemade extension
import efficientdcdist.dctree as dcdist_pascal

import time
from benchmark import create_dataset


def get_cdists(points, min_pts):
        '''
        Computes the core distances of a set of points, given a min_pts.
        '''
        num_points = points.shape[0]
        dim = int(points.shape[1])

        D = np.zeros([num_points, num_points])
        D = get_dist_matrix(points, D, dim, num_points)

        cdists = np.sort(D, axis=1)
        cdists = cdists[:, min_pts - 1] #These are the core-distances for each point.
        return cdists

@numba.njit(fastmath=True, parallel=False)
def get_dist_matrix(points, D, dim, num_points):
    '''
    Returns the Euclidean distance matrix of a 2D set of points. 

    Parameters
    ----------

    points : n x m 2D numpy array
    D : empty n x n numpy array
    dim : m
    num_points : n
    '''
    for i in numba.prange(num_points):
        x = points[i]
        for j in range(i+1, num_points):
            y = points[j]
            dist = 0
            for d in range(dim):
                dist += (x[d] - y[d]) ** 2
            dist = np.sqrt(dist)
            D[i, j] = dist
            D[j, i] = dist
    return D

k = 3
dim = 2000
n = 4000
dataset = "blobs"
points1, _ = create_dataset(num_points=n, datatype=dataset, num_features=dim)
points1t = np.transpose(points1)


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
                       [19,18]], dtype=np.float64) #If points are given as ints, we need to explicitly convert them to float.
points2t = np.transpose(points2)


points3 = np.array([[0.0, 1.0],
                    [0.0, 2.0],
                    [0.0, 3.0],
                    [0.0, 1.5],
                    [0.0, 5.5],
                    [0.0, 10.0],
                    ])
points3t = np.transpose(points3)


#dc_constructor = dcdist_pascal.DCTree(points1, k)


print("Testing dc_dist code with dataset", dataset, "and n =", n)

t1 = time.time()

cdists = dctree.compute_cdists(points1t, k, "kd")

t2 = time.time()

cdists2 = dctree.compute_cdists(points1t, k, "naive")

t3 = time.time()


#cdists3 = get_cdists(points1, k)

t4 = time.time()

cdists4 = dctree.compute_cdists(points1t, k, "naive2")

t5 = time.time()

cdists5 = dctree.compute_cdists(points1t, k, "naive3")

t6 = time.time()



print("kdtree:", t2-t1)
print("naive:", t3-t2)
print("python old", t4-t3)
print("naive2:", t5-t4)
print("naive3:", t6-t5)

#print(cdists4 - cdists3)
#print("cdist diff:", cdists-cdists3)