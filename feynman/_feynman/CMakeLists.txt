cmake_minimum_required(VERSION 3.12)

project(_feynman)

enable_language(Fortran)

# Enforce CMake to use the same version of Python, as used by setup.py script
if ("x${Python3_VERSION}" STREQUAL "x")
find_package(Python3 COMPONENTS Interpreter Development)
else()
find_package(Python3 ${Python3_VERSION} EXACT COMPONENTS Interpreter Development)
endif()
set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})

set(PYBIND11_PYTHON_VERSION ${Python3_VERSION})
add_subdirectory(ThirdParty/pybind11 EXCLUDE_FROM_ALL)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -g")

find_package(OpenMP REQUIRED)

file(GLOB SRCS "src/*.f90" "src/*.cpp")

pybind11_add_module(${PROJECT_NAME} ${SRCS})
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_Fortran)

install(TARGETS ${PROJECT_NAME} DESTINATION .)

