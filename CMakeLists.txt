cmake_minimum_required(VERSION 3.28)
include(FetchContent)
# General config
project(REEF3D LANGUAGES CXX VERSION 23.11)
# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)

# Definitions
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

FetchContent_Declare(
  openmpi
  GIT_REPOSITORY https://github.com/open-mpi/ompi.git
  GIT_TAG        v5.0.0
  CONFIGURE_COMMAND ./configure
  #FIND_PACKAGE_ARGS NAMES MPI
)
FetchContent_Declare(
    hypre
    GIT_REPOSITORY https://github.com/hypre-space/hypre.git
    GIT_TAG        v2.30.0
    SOURCE_SUBDIR src
    BINARY_DIR src/cmbuild
    CMAKE_ARGS --prefix=/usr/local/hypre
)

# MPI
find_package(MPI QUIET)
if(NOT ${MPI_FOUND})
    FetchContent_MakeAvailable(openmpi)
endif()
set(MPIEXEC_EXECUTABLE "/usr/local/openmpi/bin/mpicxx")
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})
# Set the C++ compiler to mpicxx
find_program(MPI_CXX_COMPILER mpicxx)
if(NOT MPI_CXX_COMPILER)
    message(FATAL_ERROR "mpicxx compiler not found.")
endif()

# hypre
set(HYPRE_DIR /usr/local/hypre)
if (NOT EXISTS ${HYPRE_DIR}/lib/libHYPRE.a)
    message(FATAL " hypre not found - installing")
    FetchContent_MakeAvailable(hypre)
    set(cmake_arg "../")
    set(cmake_arg2 "--install")
    set(cmake_arg3 ".")
    set(cmake_arg4 "--prefix=/usr/local/hypre")
    set(hypre_BUILD_dir /build/_deps/hypre-src/src/cmbuild)
    execute_process(COMMAND ${CMAKE_COMMAND} ${cmake_arg}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/${hypre_BUILD_dir}
    RESULT_VARIABLE hypre_result
    OUTPUT_VARIABLE hypre_ver)
    if (NOT ${hypre_ver})
        execute_process(COMMAND ${CMAKE_COMMAND} ${cmake_arg2} ${cmake_arg3}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/${hypre_BUILD_dir}
        RESULT_VARIABLE hypre_result2
        OUTPUT_VARIABLE hypre_ver2)
        message(STATUS "${hypre_ver2}")
    endif()
else()
    message(STATUS "Found hypre.")
endif()

# Eigen
set(Eigen_DIR ThirdParty/eigen-3.3.8)
set(Eigen3_DIR ${Eigen_DIR})
find_package(Eigen3 QUIET)
if(NOT ${Eigen3_FOUND})
    # Create a build directory for Eigen
    set(Eigen3_BUILD_DIR ${Eigen_DIR}/build)
    file(MAKE_DIRECTORY ${Eigen3_BUILD_DIR})
    set (cmake_arg "../")
    # Run cmake for Eigen in the build directory
    execute_process(COMMAND ${CMAKE_COMMAND} ${cmake_arg}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/${Eigen3_BUILD_DIR}
    RESULT_VARIABLE eigen_result
    OUTPUT_VARIABLE eigen_ver)
    message(STATUS "${eigen_ver}")
    set(Eigen3_DIR ${Eigen3_BUILD_DIR})
    message(STATUS "Build of Eigen3 completed.")
else()
    message(STATUS "Found Eigen3.")
endif()
find_package(Eigen3 REQUIRED)






# Compiler flags
add_library(REEF3D_compiler_flags INTERFACE)
set(gcc_like_cxx "$<COMPILE_LANG_AND_ID:CXX,ARMClang,AppleClang,Clang,GNU,LCC>")
set(msvc_cxx "$<COMPILE_LANG_AND_ID:CXX,MSVC>")
target_compile_options(REEF3D_compiler_flags INTERFACE
  "$<${gcc_like_cxx}:-w>"
  "$<${msvc_cxx}:-W3>"
)

# Globbing source files in src directory
file(GLOB_RECURSE SOURCE_FILES "src/*.cpp")

# Add executable target
add_executable(REEF3D ${SOURCE_FILES})

HEADER_DIRECTORIES(header_dir_list)

list(LENGTH header_dir_list header_dir_list_count)
message(STATUS "[INFO] Found ${header_dir_list_count} header directories.")

# Include directories
target_include_directories(REEF3D PRIVATE
    ${header_dir_list}    # Add include directories
    ${HYPRE_DIR}/include  # Add HYPRE include directory
    ${EIGEN_DIR}          # Add Eigen include directory
)

# Link directories
link_directories(
    ${HYPRE_DIR}/lib  # Add HYPRE lib directory
)

# Link libraries
target_link_libraries(REEF3D
    ${MPI_CXX_LIBRARIES}
    ${HYPRE_DIR}/lib/libHYPRE.a
    Eigen3::Eigen
    REEF3D_compiler_flags
)

# Specify output directory for the executable
set_target_properties(REEF3D PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)
