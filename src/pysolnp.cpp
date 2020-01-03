#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>

#include <vector>
#include <functional>

#include "solver.cpp"


namespace py = pybind11;

//PYBIND11_MAKE_OPAQUE(std::vector<double>);

PYBIND11_MODULE(pysolnp, m
) {
    py::class_<Result>(m, "Result")
            .def(py::init<double, std::vector<double>, cppsolnp::log_list, int>());
//    double solve_value;
//    dlib::matrix<double, 0, 1> optimum;
//    cppsolnp::log_list log;
//    int callbacks;

    m.doc() = R"pbdoc(

            -----------------------
            .. currentmodule:: pysolnp
            .. autosummary::
               :toctree: _generate
               solve
        )pbdoc";

    m.def("solve", &solve, R"pbdoc(
            Solve using the SOLNP algorithm.
        )pbdoc",
    py::arg("par_start_value"),
    py::arg("par_lower_limit"),
    py::arg("par_upper_limit"),
    py::arg("obj_func"),
    py::arg("eq_func"),
    py::arg("eq_values"),
    py::arg("ineq_func"),
    py::arg("ineq_lower_bounds"),
    py::arg("ineq_upper_bounds"),
    py::arg("rho"),
    py::arg("max_major_iter"),
    py::arg("max_minor_iter"),
    py::arg("delta"),
    py::arg("tolerance"),
    py::arg("debug"));

    #ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
    #else
    m.attr("__version__") = "dev";
    #endif
}
