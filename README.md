# Current plan: 
The parallel HBDBSCAN implementation does only work for 2-20 dimensions andhe dimension is "hard coded" in c++ templates. 
The code is also likely not efficient in higher dimensions, they only compared themselves on datasets with up to 16 dimensions.

Instead, we decided to "merge" the EMST C++ implementation with the HDBSCAN implementation in sklearn (McInnes version). 

The implementation will be done in C/C++ with python bindings.
The benchmarking will likely be done in Python.



# -----------------------------------
# Implementation steps:
1. Construct dc-tree:
    a. Core-distances
    b. MST construction
    c. The dc-tree itself from MST
2. Construct algorithms over the dc-tree:
    a. k-centroids (k-means/k-median etc.) greedy efficient algorithm
        i. New hierarchy construction
    b. HDBSCAN* over the dc-tree
    c. ???


# -----------------------------------
# Useful links for any implementation:
my thesis code: https://github.com/rsmj1/densitycenter
emst paper: https://www.mlpack.org/papers/emst.pdf
parallel hdbscan paper: https://people.csail.mit.edu/jshun/emst-hdbscan.pdf
claimed fast hdbscan (in C): https://github.com/Karthick47v2/efficient-hdbscan
- Uses quickselect to find the k-th nearest neighbor for cdists https://www.geeksforgeeks.org/quickselect-algorithm/
- Is seemingly faster on 100 dimensions but not 10. They completely switch roles on low dimensions.
- It is basically invariant of the dimensions in its speed while sklearn is not - this slows down quite a lot with higher dimensions. 

python fast hdbscan : https://github.com/TutteInstitute/fast_hdbscan?tab=readme-ov-file
state of the art fast approximate knn : https://arxiv.org/pdf/1609.07228
https://www.vldb.org/pvldb/vol12/p461-fu.pdf
https://arxiv.org/pdf/1807.05614

efanna: https://github.com/ZJULearning/efanna

state of the art of exact knn : https://www.sciencedirect.com/science/article/pii/S0950705119304678
cover tree knn (has no benchmarks) : https://papers.nips.cc/paper_files/paper/2009/hash/2421fcb1263b9530df88f7f002e78ea5-Abstract.html
dual tree concept : https://arxiv.org/pdf/1304.4327

# -----------------------------------
# Dataset loading
https://clustpy.readthedocs.io/en/latest/clustpy.data.html

These datasets are seemingly not in the above library:

KDDCUP04BIO found here: https://cs.joensuu.fi/sipu/datasets/
Taxi dataset (just the start locations): https://figshare.com/articles/dataset/Porto_taxi_trajectories/12302165
Covertype: https://archive.ics.uci.edu/dataset/31/covertype
nuScenes scene: tutorial here shows how to get a small version of it. https://www.nuscenes.org/nuscenes?tutorial=lidarseg_panoptic