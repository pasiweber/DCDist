#ifndef GRAPH_HDB_C_HPP
#define GRAPH_HDB_C_HPP

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

// gcc -shared -o hdb.so hdb_c.c -O3 -march=native -lm -fopenmp -fPIC
/*
This code computes the core distances using quickselect, 
and then computes the mutual reachability matrix from this.
*/

int compare_doubles(const void *a, const void *b);

double euclidean_distance(double *x, double *y, int dim);

void calc_distance_matrix(double *data, unsigned long long n, int dim, double *distance_matrix);

void swap(double *const a, double *const b);

unsigned long long partition(double *const arr, const unsigned long long low, const unsigned long long high);

double quickSelect(double *const arr, const unsigned long long low, const unsigned long long high, const int k);

void calc_core_dist(unsigned long long n, int k, double *core_dist, double *distance_matrix);

double *calc_mutual_reachability_dist(double *data, unsigned long long n, int dim, int k);


#endif