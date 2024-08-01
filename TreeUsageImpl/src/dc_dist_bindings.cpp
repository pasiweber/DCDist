#include <dc_dist.hpp>
#include <iostream>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
//mlpack stuff
#include <armadillo>


#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include <../parallel_hdbscan/include/hdbscan/point.h>
#include <../parallel_hdbscan/include/hdbscan/armapoint.h>
#include <../parallel_hdbscan/include/hdbscan/hdbscan.h>
#include <../parallel_hdbscan/include/hdbscan/edge.h>

#include <../parallel_hdbscan/src/kdTree.h>
#include <../parallel_hdbscan/src/kdTreeKnn.h>
#include <../parallel_hdbscan/src/kdTreeArma.h>
#include <../parallel_hdbscan/src/kdTreeKnnArma.h>
namespace py = pybind11;


using namespace parlay;


parlay::sequence<pargeo::point2> translate_points(py::array_t<double> input_array) {
    std::cout << "compute cdists wrapper called" << std::endl;
    // Convert numpy array to armadillo matrix
    auto buf = input_array.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Number of dimensions must be two");
    }

    double *ptr = static_cast<double *>(buf.ptr);
    arma::mat data(ptr, buf.shape[0], buf.shape[1], false, true);
    std::cout << "here" << std::endl;
    // Call the original compute_cdists function
    

    return convertArmaMatToParlayPoints2(data);
}



template<class T, class Seq>
py::array_t<T> wrapArray2d(Seq& result_vec, ssize_t cols) {
  ssize_t rows = (ssize_t) result_vec.size() *
    (ssize_t) sizeof(typename Seq::value_type) /
    (ssize_t) sizeof(T);
  rows /= cols;
  std::vector<ssize_t> shape = { rows, cols };
  std::vector<ssize_t> strides = { (ssize_t) sizeof(T) * cols, (ssize_t) sizeof(T) };

  return py::array(py::buffer_info(result_vec.data(),                       /* data as contiguous array  */
				   sizeof(T),                               /* size of one scalar        */
				   py::format_descriptor<T>::format(),      /* data type                 */
				   2,                                       /* number of dimensions      */
				   shape,                                   /* shape of the matrix       */
				   strides                                  /* strides for each axis     */
				   ));
}

py::array_t<double> py_hdbscan(py::array_t<double, py::array::c_style | py::array::forcecast> array, size_t minPts, std::string mode) {
  if (array.ndim() != 2)
    throw std::runtime_error("Input should be 2-D NumPy array");

    if(mode == "pargeo"){ 
        if (sizeof(pargeo::point<2>) != 16)
            throw std::runtime_error("sizeof(pargeo::point<2>) != 16, check point.h");

        int dim = array.shape()[1];
        std::cout << "dim:" << dim << std::endl;
        size_t n = array.size() / dim;
        parlay::sequence<pargeo::dirEdge> result_vec;

        sequence<pargeo::wghEdge> E;
        if (dim == 2) {
            parlay::sequence<pargeo::point<2>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<2>(P, minPts);
        } else if (dim == 20) {
            parlay::sequence<pargeo::point<20>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<20>(P, minPts);
        } else {
            throw std::runtime_error("Only dimensions 2-20 is supported at the moment");
        }

        sequence<pargeo::dendroNode> dendro = pargeo::dendrogram(E, n);
        sequence<double> A(dendro.size()*4);
        parlay::parallel_for(0, dendro.size(), [&](size_t i){
                            A[i*4+0] = std::get<0>(dendro[i]);
                            A[i*4+1] = std::get<1>(dendro[i]);
                            A[i*4+2] = std::get<2>(dendro[i]);
                            A[i*4+3] = std::get<3>(dendro[i]);});
        return wrapArray2d<double>(A, 4);
    } else{
        //TODO: This seemingly does not work for very few points?? Like 10 points
        parlay::sequence<pargeo::point2> points = translate_points(array);
        int n = points.size();
        int dim = points[0].dim; 
        std::cout << "n: " << n << ", dim: " << dim << std::endl;
        parlay::sequence<pargeo::dirEdge> result_vec;
        sequence<pargeo::wghEdge> E;
        E = pargeo::hdbscan_arma(points, minPts);

        sequence<pargeo::dendroNode> dendro = pargeo::dendrogram(E, n);
        sequence<double> A(dendro.size()*4);
        parlay::parallel_for(0, dendro.size(), [&](size_t i){
                            A[i*4+0] = std::get<0>(dendro[i]);
                            A[i*4+1] = std::get<1>(dendro[i]);
                            A[i*4+2] = std::get<2>(dendro[i]);
                            A[i*4+3] = std::get<3>(dendro[i]);});

        return wrapArray2d<double>(A, 4);
    }
}




std::vector<double> compute_cdists_wrapper(py::array_t<double> input_array, size_t k, std::string mode) {
    std::cout << "compute cdists wrapper called" << std::endl;
    // Convert numpy array to armadillo matrix
    auto buf = input_array.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Number of dimensions must be two");
    }

    double *ptr = static_cast<double *>(buf.ptr);
    arma::mat data(ptr, buf.shape[0], buf.shape[1], false, true);
    std::cout << "here" << std::endl;
    // Call the original compute_cdists function
    std::vector<double> result = compute_cdists(data, k, mode);

    return result;
}

PYBIND11_MODULE(dctree, m) {
    m.doc() = "Bindings for compute_cdists function"; // Optional module docstring

    m.def("compute_cdists", &compute_cdists_wrapper, "Compute cdists",
          py::arg("data"), py::arg("k"), py::arg("mode"));

    m.def("compute_hdbscan_labels", &py_hdbscan, "Compute HDBSCAN and return integers",
          py::arg("array"), py::arg("minPts"), py::arg("mode"));
}
