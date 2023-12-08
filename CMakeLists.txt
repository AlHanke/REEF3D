cmake_minimum_required(VERSION 3.10)
project(REEF3D LANGUAGES CXX VERSION 23.11)

MACRO(HEADER_DIRECTORIES return_list)
    FILE(GLOB_RECURSE new_list src/*.h)
    SET(dir_list "")
    FOREACH(file_path ${new_list})
        GET_FILENAME_COMPONENT(dir_path ${file_path} PATH)
        SET(dir_list ${dir_list} ${dir_path})
    ENDFOREACH()
    LIST(REMOVE_DUPLICATES dir_list)
    SET(${return_list} ${dir_list})
ENDMACRO()

set(MPIEXEC_EXECUTABLE "/usr/local/openmpi/bin/mpicxx")
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
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
set(Eigen3_DIR ${Eigen3_BUILD_DIR})
find_package(Eigen3 3.3 REQUIRED)

# Globbing source files in src directory
file(GLOB_RECURSE SOURCE_FILES "src/*.cpp")

# Add executable target
add_executable(REEF3D ${SOURCE_FILES})

HEADER_DIRECTORIES(header_dir_list)

list(LENGTH header_dir_list header_dir_list_count)
message(STATUS "[INFO] Found ${header_dir_list_count} header directories.")

# Include directories
target_include_directories(REEF3D PRIVATE
    include
    ${header_dir_list} # Recursive
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
