cmake_minimum_required(VERSION 3.1)

project(_feynman)

enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")

add_executable(${PROJECT_NAME}_symbolic_regress1 symbolic_regress1.f90)
add_executable(${PROJECT_NAME}_symbolic_regress2 symbolic_regress2.f90)
add_executable(${PROJECT_NAME}_symbolic_regress3 symbolic_regress3.f90)

install(TARGETS ${PROJECT_NAME}_symbolic_regress1 DESTINATION .)
install(TARGETS ${PROJECT_NAME}_symbolic_regress2 DESTINATION .)
install(TARGETS ${PROJECT_NAME}_symbolic_regress3 DESTINATION .)
