#include <pybind11/pybind11.h>

extern "C" void symbolic_regress1();
extern "C" void symbolic_regress2();
extern "C" void symbolic_regress3();

PYBIND11_MODULE(_feynman, m) {
    m.doc() = "Feynman symbolic regression native module";
    m.def("symbolic_regress1", &symbolic_regress1);
    m.def("symbolic_regress2", &symbolic_regress2);
    m.def("symbolic_regress3", &symbolic_regress3);
}

