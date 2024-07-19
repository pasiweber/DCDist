#include <dc_dist.hpp>
#include <iostream>
#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
//mlpack stuff
#include <mlpack.hpp>



namespace py = pybind11;



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
}
