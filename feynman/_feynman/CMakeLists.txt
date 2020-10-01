cmake_minimum_required(VERSION 3.12)

project(_feynman)

enable_language(Fortran)

#set(CMAKE_BUILD_TYPE Debug)

if (UNIX)
# Use $ORIGIN or @loader_path, which allow to pick up
# dependencies relative to the executable folder
if (APPLE)
set(CMAKE_INSTALL_RPATH "@loader_path")
else()
set(CMAKE_INSTALL_RPATH "$ORIGIN")
endif()
endif()

# Enforce CMake to use the same version of Python, as used by setup.py script
if ("x${Python3_VERSION}" STREQUAL "x")
find_package(Python3 COMPONENTS Interpreter Development)
else()
find_package(Python3 ${Python3_VERSION} EXACT COMPONENTS Interpreter Development)
endif()
set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE})

set(PYBIND11_PYTHON_VERSION ${Python3_VERSION})
add_subdirectory(ThirdParty/pybind11 EXCLUDE_FROM_ALL)

set(TBB_BUILD_TESTS OFF CACHE BOOL "" FORCE)
add_subdirectory(ThirdParty/tbb EXCLUDE_FROM_ALL)
set_target_properties(tbb PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(tbbmalloc PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set_target_properties(tbbmalloc_proxy PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
set(TBB_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/tbb/include)
set(TBB_LIBRARIES tbb tbbmalloc tbbmalloc_proxy)
add_library(Thrust::TBB INTERFACE IMPORTED)
set_target_properties(Thrust::TBB PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${TBB_INCLUDE_DIR}")
target_link_libraries(Thrust::TBB INTERFACE "${TBB_LIBRARIES}")

set(THRUST_REQUIRED_SYSTEMS TBB)
add_subdirectory(ThirdParty/thrust EXCLUDE_FROM_ALL)
_thrust_setup_system(TBB)
thrust_create_target(ThrustTBB HOST TBB DEVICE TBB)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -g")

file(GLOB SRCS "src/*.f90" "src/*.cpp")

pybind11_add_module(${PROJECT_NAME} ${SRCS})
target_link_libraries(${PROJECT_NAME} PRIVATE Thrust::Thrust ThrustTBB)

install(TARGETS ${PROJECT_NAME} DESTINATION .)
install(FILES $<TARGET_FILE:tbb> DESTINATION .)
install(FILES $<TARGET_FILE:tbbmalloc> DESTINATION .)
install(FILES $<TARGET_FILE:tbbmalloc_proxy> DESTINATION .)
