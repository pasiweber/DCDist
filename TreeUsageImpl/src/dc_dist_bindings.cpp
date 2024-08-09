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



std::pair<double*, std::pair<size_t, size_t>> load_data(py::array_t<double> array){
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
    std::cout << "dim: " << cols << ", n: " << rows << std::endl;

    return {data, {rows, cols}};

}


std::vector<std::vector<int>> get_c_cluster_solutions(py::array_t<double> array, size_t minPts, size_t mcs, size_t k, py::list modes) {
    /*
        1. Create the tree
        2. Provided list of string mode inputs that in a big if else in a for loop executes the different things over the tree
            ,  resulting in different labels that will be appended and returned
    */
    std::pair<double*, std::pair<size_t, size_t>> array_data = load_data(array);
    double* data = array_data.first;
    int dim = array_data.second.second;
    unsigned long long n = array_data.second.first;
    // std::cout << "dim: " << dim << ", n: " << n << std::endl;

    // for(int i = 0; i < n*dim; i++){
    //     std::cout << " " << data[i];
    // }

    Node* dc_tree = construct_dc_tree(data, n, dim, minPts);

    std::vector<std::vector<int>> result_labels;
    result_labels.push_back({0,1,2});
    result_labels.push_back({9,8,7});

    for(const auto& item : modes){
        std::string mode = item.cast<std::string>();
        std::cout << "mode: " << mode << std::endl;
    }


    return result_labels;
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


    //py::array_t<double> array, size_t minPts, size_t mcs, size_t k, py::array_t<std::string> modes
    m.def("compute_clustering", &get_c_cluster_solutions, "Compute HDBSCAN and return integers",
          py::arg("array"), py::arg("minPts"), py::arg("mcs"), py::arg("k"), py::arg("modes"));
}
