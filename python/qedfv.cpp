#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <qedfv/coef.hpp>

namespace py = pybind11;
using namespace qedfv;
using namespace pybind11::literals;

PYBIND11_MODULE(pyqedfv, m)
{
  py::class_<Coef::Params>(m, "CoefParameters")
      .def(py::init<>())
      .def_readwrite("eta", &Coef::Params::eta)
      .def_readwrite("nmax", &Coef::Params::nmax)
      .def("__repr__",
           [](const Coef::Params &p) {
             return "{ eta: " + std::to_string(p.eta) + ", nmax: " + std::to_string(p.nmax) + " }";
           });

  py::enum_<Coef::Qed>(m, "Qed").value("L", Coef::Qed::L).value("r", Coef::Qed::r);

  py::class_<Coef>(m, "Coef")
      .def(py::init<const Coef::Qed, const bool>(), "qed"_a = Coef::Qed::L,
           "debug"_a = false)
      .def("__call__",
           static_cast<double (Coef::*)(const double, const double, const unsigned int)>(
               &Coef::operator()))
      .def("__call__",
           static_cast<double (Coef::*)(const double, const DVec3, const double,
                                             const unsigned int)>(&Coef::operator()))
      .def("__call__", static_cast<double (Coef::*)(const double, const Coef::Params)>(
                           &Coef::operator()))
      .def("__call__",
           static_cast<double (Coef::*)(const double, const DVec3, const Coef::Params)>(
               &Coef::operator()))
      .def("a", &Coef::a)
      .def("r", &Coef::r)
      .def("rBar", &Coef::rBar)
      .def("tune",
           static_cast<Coef::Params (Coef::*)(const double, const double, const double,
                                                        const double, const unsigned int,
                                                        const unsigned int)>(&Coef::tune),
           "j"_a, "residual"_a = QEDFV_DEFAULT_ERROR, "eta0"_a = 1.0, "etaFactor"_a = 0.98,
           "nmax0"_a = 5, "nmaxStep"_a = 5)
      .def("tune",
           static_cast<Coef::Params (Coef::*)(const double, const Coef::Params,
                                                        const double, const double,
                                                        const unsigned int)>(&Coef::tune),
           "j"_a, "par"_a, "residual"_a = QEDFV_DEFAULT_ERROR, "etaFactor"_a = 0.98,
           "nmaxStep"_a = 5)
      .def("tune",
           static_cast<Coef::Params (Coef::*)(
               const double, const DVec3, const double, const double, const double,
               const unsigned int, const unsigned int)>(&Coef::tune),
           "j"_a, "v"_a, "residual"_a = QEDFV_DEFAULT_ERROR, "eta0"_a = 1.0, "etaFactor"_a = 0.98,
           "nmax0"_a = 5, "nmaxStep"_a = 5)
      .def("tune",
           static_cast<Coef::Params (Coef::*)(
               const double, const DVec3, const Coef::Params, const double, const double,
               const unsigned int)>(&Coef::tune),
           "j"_a, "v"_a, "par"_a, "residual"_a = QEDFV_DEFAULT_ERROR, "etaFactor"_a = 0.98,
           "nmaxStep"_a = 5);
}