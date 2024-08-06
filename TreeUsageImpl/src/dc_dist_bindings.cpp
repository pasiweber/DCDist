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
#include <../parallel_hdbscan/include/hdbscan/vectorpoint.h>

#include <../parallel_hdbscan/include/hdbscan/hdbscan.h>
#include <../parallel_hdbscan/include/hdbscan/edge.h>

#include <../parallel_hdbscan/src/kdTree.h>
#include <../parallel_hdbscan/src/kdTreeKnn.h>
#include <../parallel_hdbscan/src/kdTreeArma.h>
#include <../parallel_hdbscan/src/kdTreeKnnArma.h>
namespace py = pybind11;


using namespace parlay;


parlay::sequence<pargeo::ArmaPoint> convertArmaMatToParlayPoints_Arma(std::vector<std::vector<double>> mat) {
    // Create a parlay sequence of points
    parlay::sequence<pargeo::ArmaPoint> points(mat.size());
    std::cout<< "num points:" << mat.size() << std::endl;
    // Populate the sequence with points created from the matrix columns
    for (size_t i = 0; i < mat.size(); ++i) {
        arma::Col<double> coords(mat[0].size());
        for (size_t j = 0; j < mat[0].size(); ++j) {
            coords[j] = mat[i][j];
        }
        points[i] = pargeo::ArmaPoint(&coords);
    }
    return points;
}

parlay::sequence<pargeo::VectorPoint> convertArmaMatToParlayPoints_Vector(std::vector<std::vector<double>> mat) {
    // Create a parlay sequence of points
    parlay::sequence<pargeo::VectorPoint> points(mat.size());
    std::cout<< "num points:" << mat.size() << std::endl;
    // Populate the sequence with points created from the matrix columns
    for (size_t i = 0; i < mat.size(); ++i) {
        std::vector<double> coords(mat[0].size());
        for (size_t j = 0; j < mat[0].size(); ++j) {
            coords[j] = mat[i][j];
        }
        points[i] = pargeo::VectorPoint(&coords);
    }
    return points;
}
parlay::sequence<pargeo::ArmaPoint> translate_points_arma(py::array_t<double> input_array) {
    std::cout << "compute cdists wrapper called" << std::endl;
    // Convert numpy array to armadillo matrix
    auto buf = input_array.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Number of dimensions must be two");
    }

    double *ptr = static_cast<double *>(buf.ptr);
    ssize_t rows = buf.shape[0];
    ssize_t cols = buf.shape[1];
    std::cout << "rows, cols:" << rows << ", " << cols << std::endl;

    std::vector<std::vector<double>> result(rows, std::vector<double>(cols));

    // Populate the vector of vectors
    for (ssize_t i = 0; i < rows; ++i) {
        for (ssize_t j = 0; j < cols; ++j) {
            result[i][j] = ptr[i * cols + j];
        }
    }
    std::cout << "here" << std::endl;
    // Call the original compute_cdists function
    return convertArmaMatToParlayPoints_Arma(result);
}

parlay::sequence<pargeo::VectorPoint> translate_points_vector(py::array_t<double> input_array) {
    std::cout << "compute cdists wrapper called" << std::endl;
    // Convert numpy array to armadillo matrix
    auto buf = input_array.request();
    if (buf.ndim != 2) {
        throw std::runtime_error("Number of dimensions must be two");
    }

    double *ptr = static_cast<double *>(buf.ptr);
    ssize_t rows = buf.shape[0];
    ssize_t cols = buf.shape[1];
    std::cout << "rows, cols:" << rows << ", " << cols << std::endl;

    std::vector<std::vector<double>> result(rows, std::vector<double>(cols));

    // Populate the vector of vectors
    for (ssize_t i = 0; i < rows; ++i) {
        for (ssize_t j = 0; j < cols; ++j) {
            result[i][j] = ptr[i * cols + j];
        }
    }
    std::cout << "here" << std::endl;
    // Call the original compute_cdists function
    return convertArmaMatToParlayPoints_Vector(result);
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
            // for(int i = 0; i < P.size(); i++){
            //     std::cout << "i: " << i << ", p: " << P[i] << std::endl;
            // }

            E = pargeo::hdbscan<2>(P, minPts);
        } else if (dim == 200) {
            parlay::sequence<pargeo::point<200>> P(n);
            std::memcpy(P.data(), array.data(), array.size() * sizeof(double));

            // for(int i = 0; i < P.size(); i++){
            //     std::cout << "i: " << i << ", p: " << P[i] << std::endl;
            // }


            E = pargeo::hdbscan<200>(P, minPts);
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
    } else if(mode == "arma"){
        std::cout << "Arma version" << std::endl;
        //TODO: This seemingly does not work for very few points?? Like 10 points
        parlay::sequence<pargeo::ArmaPoint> points = translate_points_arma(array);

        // for(int i = 0; i < points.size(); i++){
        //     std::cout << "i: " << i << ", p: " << points[i] << std::endl;
        // }

        


        int n = points.size();
        int dim = points[0].dim; 
        std::cout << "n: " << n << ", dim: " << dim << std::endl;
        parlay::sequence<pargeo::dirEdge> result_vec;
        sequence<pargeo::wghEdge> E;
        E = pargeo::hdbscan_arma(points, minPts);
        
        std::cout << "Done with HDBSCAN part, Computing dendrogram" << std::endl;
        sequence<pargeo::dendroNode> dendro = pargeo::dendrogram(E, n);
        sequence<double> A(dendro.size()*4);
        parlay::parallel_for(0, dendro.size(), [&](size_t i){
                            A[i*4+0] = std::get<0>(dendro[i]);
                            A[i*4+1] = std::get<1>(dendro[i]);
                            A[i*4+2] = std::get<2>(dendro[i]);
                            A[i*4+3] = std::get<3>(dendro[i]);});

        return wrapArray2d<double>(A, 4);
    } else {
        std::cout << "Vector version" << std::endl;
        //TODO: This seemingly does not work for very few points?? Like 10 points
        parlay::sequence<pargeo::VectorPoint> points = translate_points_vector(array);
        int n = points.size();
        int dim = points[0].dim; 
        std::cout << "n: " << n << ", dim: " << dim << std::endl;
        parlay::sequence<pargeo::dirEdge> result_vec;
        sequence<pargeo::wghEdge> E;
        E = pargeo::hdbscan_vector(points, minPts);

        std::cout << "Done with HDBSCAN part, Computing dendrogram" << std::endl;
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
