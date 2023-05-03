#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <qedfvcoef.hpp>

namespace py = pybind11;
using namespace qedfv;
using namespace pybind11::literals;

PYBIND11_MODULE(qedfv, m)
{
  py::class_<QedFvCoef::Params>(m, "CoefParameters")
      .def(py::init<>())
      .def_readwrite("eta", &QedFvCoef::Params::eta)
      .def_readwrite("nmax", &QedFvCoef::Params::nmax)
      .def("__repr__",
           [](const QedFvCoef::Params &p) {
             return "{ eta: " + std::to_string(p.eta) + ", nmax: " + std::to_string(p.nmax) + " }";
           });

  py::enum_<QedFvCoef::Qed>(m, "Qed").value("L", QedFvCoef::Qed::L).value("r", QedFvCoef::Qed::r);

  py::class_<QedFvCoef>(m, "Coef")
      .def(py::init<const QedFvCoef::Qed, const bool>(), "qed"_a = QedFvCoef::Qed::L,
           "debug"_a = false)
      .def("__call__",
           static_cast<double (QedFvCoef::*)(const double, const double, const unsigned int)>(
               &QedFvCoef::operator()))
      .def("__call__",
           static_cast<double (QedFvCoef::*)(const double, const DVec3, const double,
                                             const unsigned int)>(&QedFvCoef::operator()))
      .def("__call__", static_cast<double (QedFvCoef::*)(const double, const QedFvCoef::Params)>(
                           &QedFvCoef::operator()))
      .def("__call__",
           static_cast<double (QedFvCoef::*)(const double, const DVec3, const QedFvCoef::Params)>(
               &QedFvCoef::operator()))
      .def("a", &QedFvCoef::a)
      .def("r", &QedFvCoef::r)
      .def("rBar", &QedFvCoef::rBar)
      .def("tune",
           static_cast<QedFvCoef::Params (QedFvCoef::*)(const double, const double, const double,
                                                        const double, const unsigned int,
                                                        const unsigned int)>(&QedFvCoef::tune),
           "j"_a, "residual"_a = QEDFV_DEFAULT_ERROR, "eta0"_a = 1.0, "etaFactor"_a = 0.98,
           "nmax0"_a = 5, "nmaxStep"_a = 5)
      .def("tune",
           static_cast<QedFvCoef::Params (QedFvCoef::*)(
               const double, const DVec3, const double, const double, const double,
               const unsigned int, const unsigned int)>(&QedFvCoef::tune),
           "j"_a, "v"_a, "residual"_a = QEDFV_DEFAULT_ERROR, "eta0"_a = 1.0, "etaFactor"_a = 0.98,
           "nmax0"_a = 5, "nmaxStep"_a = 5);
}