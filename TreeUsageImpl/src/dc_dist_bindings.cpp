#include <dc_dist.hpp>
#include <iostream>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
//mlpack stuff
#include <armadillo>


#include <dc_hdbscan.hpp>

namespace py = pybind11;




// parlay::sequence<pargeo::ArmaPoint> convertArmaMatToParlayPoints_Arma(std::vector<std::vector<double>> mat) {
//     // Create a parlay sequence of points
//     parlay::sequence<pargeo::ArmaPoint> points(mat.size());
//     std::cout<< "num points:" << mat.size() << std::endl;
//     // Populate the sequence with points created from the matrix columns
//     for (size_t i = 0; i < mat.size(); ++i) {
//         arma::Col<double> coords(mat[0].size());
//         for (size_t j = 0; j < mat[0].size(); ++j) {
//             coords[j] = mat[i][j];
//         }
//         points[i] = pargeo::ArmaPoint(&coords);
//     }
//     return points;
// }

// parlay::sequence<pargeo::VectorPoint> convertArmaMatToParlayPoints_Vector(std::vector<std::vector<double>> mat) {
//     // Create a parlay sequence of points
//     parlay::sequence<pargeo::VectorPoint> points(mat.size());
//     std::cout<< "num points:" << mat.size() << std::endl;
//     // Populate the sequence with points created from the matrix columns
//     for (size_t i = 0; i < mat.size(); ++i) {
//         std::vector<double> coords(mat[0].size());
//         for (size_t j = 0; j < mat[0].size(); ++j) {
//             coords[j] = mat[i][j];
//         }
//         points[i] = pargeo::VectorPoint(&coords);
//     }
//     return points;
// }
// parlay::sequence<pargeo::ArmaPoint> translate_points_arma(py::array_t<double> input_array) {
//     std::cout << "compute cdists wrapper called" << std::endl;
//     // Convert numpy array to armadillo matrix
//     auto buf = input_array.request();
//     if (buf.ndim != 2) {
//         throw std::runtime_error("Number of dimensions must be two");
//     }

//     double *ptr = static_cast<double *>(buf.ptr);
//     ssize_t rows = buf.shape[0];
//     ssize_t cols = buf.shape[1];
//     std::cout << "rows, cols:" << rows << ", " << cols << std::endl;

//     std::vector<std::vector<double>> result(rows, std::vector<double>(cols));

//     // Populate the vector of vectors
//     for (ssize_t i = 0; i < rows; ++i) {
//         //std::cout << "i " << i << " : ";
//         for (ssize_t j = 0; j < cols; ++j) {
//             result[i][j] = ptr[i * cols + j];
//             //std::cout << result[i][j] << " ";
//         }
//         //std::cout << std::endl;

//     }
//     //std::cout << "here" << std::endl;
//     // Call the original compute_cdists function
//     return convertArmaMatToParlayPoints_Arma(result);
// }

// parlay::sequence<pargeo::VectorPoint> translate_points_vector(py::array_t<double> input_array) {
//     std::cout << "compute cdists wrapper called" << std::endl;
//     // Convert numpy array to armadillo matrix
//     auto buf = input_array.request();
//     if (buf.ndim != 2) {
//         throw std::runtime_error("Number of dimensions must be two");
//     }

//     double *ptr = static_cast<double *>(buf.ptr);
//     ssize_t rows = buf.shape[0];
//     ssize_t cols = buf.shape[1];
//     std::cout << "rows, cols:" << rows << ", " << cols << std::endl;

//     std::vector<std::vector<double>> result(rows, std::vector<double>(cols));

//     // Populate the vector of vectors
//     for (ssize_t i = 0; i < rows; ++i) {
//         for (ssize_t j = 0; j < cols; ++j) {
//             result[i][j] = ptr[i * cols + j];
//         }
//     }
//     std::cout << "here" << std::endl;
//     // Call the original compute_cdists function
//     return convertArmaMatToParlayPoints_Vector(result);
// }

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

std::vector<int> py_hdbscan(py::array_t<double> array, size_t minPts, size_t mcs) {
    py::buffer_info buf_info = array.request();

    if (buf_info.ndim != 2) {
        throw std::runtime_error("Input array must be 2-dimensional");
    }
    // Get the shape of the array
    size_t rows = buf_info.shape[0];
    size_t cols = buf_info.shape[1];
    size_t size = rows * cols;
    double* data = new double[size];

    std::memcpy(data, buf_info.ptr, size * sizeof(double));

    int dim = cols;
    unsigned long long n = rows;
    std::cout << "dim: " << cols << ", n: " << n << std::endl;

    Dc_hdbscan tree(minPts, mcs);
    tree.fit(data, n, dim, minPts);

    std::vector<int> labels = tree.labels_;
    return labels;
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
          py::arg("array"), py::arg("minPts"), py::arg("mcs"));
}
