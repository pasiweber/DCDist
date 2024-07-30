#include <dc_dist.hpp>
#include <iostream>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
//mlpack stuff
#include <mlpack.hpp>


#include "parlay/parallel.h"
#include "parlay/sequence.h"
#include <../parallel_hdbscan/src/kdTree.h>
#include <../parallel_hdbscan/src/kdTreeKnn.h>
#include <../parallel_hdbscan/include/hdbscan/point.h>
#include <../parallel_hdbscan/src/kdTreeArma.h>
#include <../parallel_hdbscan/src/kdTreeKnnArma.h>
#include <../parallel_hdbscan/include/hdbscan/armapoint.h>
#include <../parallel_hdbscan/include/hdbscan/hdbscan.h>
namespace py = pybind11;



using namespace parlay;

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
        size_t n = array.size() / dim;
        parlay::sequence<pargeo::dirEdge> result_vec;

        sequence<pargeo::wghEdge> E;
        if (dim == 2) {
            parlay::sequence<pargeo::point<2>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<2>(P, minPts);
        } else if (dim == 3) {
            parlay::sequence<pargeo::point<3>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<3>(P, minPts);
        } else if (dim == 4) {
            parlay::sequence<pargeo::point<4>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<4>(P, minPts);
        } else if (dim == 5) {
            parlay::sequence<pargeo::point<5>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<5>(P, minPts);
        } else if (dim == 6) {
            parlay::sequence<pargeo::point<6>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<6>(P, minPts);
        } else if (dim == 7) {
            parlay::sequence<pargeo::point<7>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<7>(P, minPts);
        } else if (dim == 8) {
            parlay::sequence<pargeo::point<8>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<8>(P, minPts);
        } else if (dim == 9) {
            parlay::sequence<pargeo::point<9>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<9>(P, minPts);
        } else if (dim == 10) {
            parlay::sequence<pargeo::point<10>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<10>(P, minPts);
        } else if (dim == 11) {
            parlay::sequence<pargeo::point<11>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<11>(P, minPts);
        } else if (dim == 12) {
            parlay::sequence<pargeo::point<12>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<12>(P, minPts);
        } else if (dim == 13) {
            parlay::sequence<pargeo::point<13>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<13>(P, minPts);
        } else if (dim == 14) {
            parlay::sequence<pargeo::point<14>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<14>(P, minPts);
        } else if (dim == 15) {
            parlay::sequence<pargeo::point<15>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<15>(P, minPts);
        } else if (dim == 16) {
            parlay::sequence<pargeo::point<16>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<16>(P, minPts);
        } else if (dim == 17) {
            parlay::sequence<pargeo::point<17>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<17>(P, minPts);
        } else if (dim == 18) {
            parlay::sequence<pargeo::point<18>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<18>(P, minPts);
        } else if (dim == 19) {
            parlay::sequence<pargeo::point<19>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
            E = pargeo::hdbscan<19>(P, minPts);
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

        int dim = array.shape()[1];
        size_t n = array.size() / dim;
        parlay::sequence<pargeo::dirEdge> result_vec;

        sequence<pargeo::wghEdge> E;

        parlay::sequence<pargeo::point2> P(n); //This might not work, we'll see
        std::memcpy(P.data(), array.data(), array.size() * sizeof(double));
        E = pargeo::hdbscan_arma(P, minPts);
        
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
