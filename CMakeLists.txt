cmake_minimum_required(VERSION 3.2)
project(MuEOSLibrary VERSION 1.0 LANGUAGES C CXX)
#add_definitions(-DSOURCE_DIR="${CMAKE_SOURCE_DIR}")

find_package(HDF5 REQUIRED COMPONENTS C CXX HL)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # Optional: disallow falling back to older standards
set(CMAKE_CXX_EXTENSIONS OFF)        # Optional: disable compiler-specific extensions (like gnu++17)

#set(CMAKE_CXX_FLAGS "-O3 -Did_test=0")

set(IDX_EOS 1 CACHE STRING "Index for lepton EOS mode")

include_directories(
    src
    src/eos_species/
)

################################################################################

macro(ADD_MUEOSLIBRARY_EXECUTABLE NAME PATH ID)
add_executable(
    ${NAME}
    ${PATH}
    ${PROJECT_SOURCE_DIR}/src/eos_species/eos_species.cpp
    ${PROJECT_SOURCE_DIR}/src/eos_species/eos_neutrinos.cpp
    ${PROJECT_SOURCE_DIR}/src/eos_species/eos_photons.cpp
    ${PROJECT_SOURCE_DIR}/src/eos_species/eos_baryons.cpp
    ${PROJECT_SOURCE_DIR}/src/helmholtz_eos/helmholtz_eos.cpp
    ${PROJECT_SOURCE_DIR}/src/fermi_integrals/nonrel_FD.cpp
    ${PROJECT_SOURCE_DIR}/src/fermi_integrals/generalized_FD.cpp
    ${PROJECT_SOURCE_DIR}/src/num_tools/safe_exp.cpp
    ${PROJECT_SOURCE_DIR}/src/num_tools/root_finding/newt_raphson.cpp
)
target_compile_definitions(${NAME} PUBLIC id_test=${ID})
target_include_directories(${NAME} PRIVATE ${HDF5_INCLUDE_DIRS})
target_link_libraries(${NAME}
    m
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}
)
endmacro()

################################################################################

ADD_MUEOSLIBRARY_EXECUTABLE(
main
main.cpp
${IDX_EOS})

ADD_MUEOSLIBRARY_EXECUTABLE(
output_table
output/generate_output_table.cpp
2)

#ADD_MUEOSLIBRARY_EXECUTABLE(
#generate_data_plot
#output/generate_data_plot.cpp
#2)

ADD_MUEOSLIBRARY_EXECUTABLE(
lep_table
output/generate_lep_table.cpp
2)

#ADD_MUEOSLIBRARY_EXECUTABLE(
#comparison_with_mu
#tests/src/comparison_with_mu.cpp
#2)

#ADD_MUEOSLIBRARY_EXECUTABLE(
#comparison_without_mu
#tests/src/comparison_without_mu.cpp
#2)

#ADD_MUEOSLIBRARY_EXECUTABLE(
#hist_comp_with_mu
#tests/src/hist_comp_with_mu.cpp
#2)

#ADD_MUEOSLIBRARY_EXECUTABLE(
#number_nr_iter
#tests/src/number_nr_iter.cpp
#2)
