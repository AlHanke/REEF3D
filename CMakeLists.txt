cmake_minimum_required(VERSION 3.10)
project(REEF3D CXX VERSION 23.11)

# Set the C++ compiler to mpicxx
find_program(MPI_CXX_COMPILER mpicxx)
if(NOT MPI_CXX_COMPILER)
    message(FATAL_ERROR "mpicxx compiler not found.")
endif()

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)

add_library(REEF3D_compiler_flags INTERFACE)
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
target_compile_options(REEF3D_compiler_flags INTERFACE
  "$<${gcc_like_cxx}:-w>"
  "$<${msvc_cxx}:-W3>"
)

set(HYPRE_DIR /usr/local/hypre)
set(Eigen_DIR ThirdParty/eigen-3.3.8)

# Create a build directory for Eigen
set(Eigen3_BUILD_DIR ${Eigen_DIR}/build)
file(MAKE_DIRECTORY ${Eigen3_BUILD_DIR})

# Run cmake for Eigen in the build directory
execute_process(
    COMMAND cmake ${CMAKE_SOURCE_DIR}/${Eigen3_BUILD_DIR}
    WORKING_DIRECTORY ${Eigen3_BUILD_DIR}
)
find_package(Eigen3 3.3 REQUIRED)
set(Eigen3_DIR T${Eigen3_BUILD_DIR})

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Globbing source files in src directory
file(GLOB SOURCE_FILES "src/*.cpp")

# Add executable target
add_executable(REEF3D ${SOURCE_FILES})

# Include directories
target_include_directories(REEF3D PRIVATE
    include
    ${HYPRE_DIR}/include  # Add your HYPRE include directory
    ${EIGEN_DIR}          # Add your Eigen include directory
)

# Link directories
link_directories(
    ${HYPRE_DIR}/lib  # Add your HYPRE lib directory
)

# Link libraries
target_link_libraries(REEF3D
    ${MPI_CXX_LIBRARIES}
    Eigen3::Eigen
    ${HYPRE_DIR}/lib/libHYPRE.a
    REEF3D_compiler_flags
)

# Specify output directory for the executable
set_target_properties(REEF3D PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
