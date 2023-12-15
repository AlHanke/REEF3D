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
    hypre
    GIT_REPOSITORY https://github.com/hypre-space/hypre.git
    GIT_TAG        v2.30.0
    SOURCE_SUBDIR src
    BINARY_DIR cmbuild
)


# MPI
find_package(MPI QUIET)
if(NOT ${MPI_FOUND})
    # Install openMPI
    message(FATAL " MPI not found - installing openMPI")
    set(OPENMPI_VERSION_MAJOR  "5")
    set(OPENMPI_VERSION_MINOR  "0")
    set(OPENMPI_VERSION_PATCH  "0")
    set(OPENMPI_DIR _deps/openmpi-${OPENMPI_VERSION_MAJOR}.${OPENMPI_VERSION_MINOR}.${OPENMPI_VERSION_PATCH})
    if (NOT EXISTS ${OPENMPI_DIR})
        # Download openMPI
        file(DOWNLOAD
            https://download.open-mpi.org/release/open-mpi/v${OPENMPI_VERSION_MAJOR}.${OPENMPI_VERSION_MINOR}/openmpi-${OPENMPI_VERSION_MAJOR}.${OPENMPI_VERSION_MINOR}.${OPENMPI_VERSION_PATCH}.tar.bz2
            ${PROJECT_BINARY_DIR}/_deps/openmpi-${OPENMPI_VERSION_MAJOR}.${OPENMPI_VERSION_MINOR}.${OPENMPI_VERSION_PATCH}.tar.bz2
            EXPECTED_MD5 5c4fae19e1d853abeb1bdadd10a59466
        )
        execute_process(COMMAND ${CMAKE_COMMAND} -E tar xf openmpi-${OPENMPI_VERSION_MAJOR}.${OPENMPI_VERSION_MINOR}.${OPENMPI_VERSION_PATCH}.tar.bz2
                        WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/_deps
                        RESULT_VARIABLE ompi_result
                        OUTPUT_VARIABLE ompi_ver)
        if (${ompi_result}==1)
            message("${ompi_ver}")
        endif()
        file(REMOVE ${PROJECT_BINARY_DIR}/_deps/openmpi-${OPENMPI_VERSION_MAJOR}.${OPENMPI_VERSION_MINOR}.${OPENMPI_VERSION_PATCH}.tar.bz2)
    endif()
    execute_process(COMMAND ./configure --prefix=/usr/local/openmpi
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${OPENMPI_DIR}
                    RESULT_VARIABLE ompi_result
                    OUTPUT_VARIABLE ompi_ver)
    if (${ompi_result}==1)
        message("${ompi_ver}")
    endif()
    cmake_host_system_information(RESULT N
    QUERY NUMBER_OF_PHYSICAL_CORES)
    execute_process(COMMAND make -j ${N} all -w
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${OPENMPI_DIR}
                    RESULT_VARIABLE ompi_result
                    OUTPUT_VARIABLE ompi_ver)
    if (${ompi_result}==1)
        message("${ompi_ver}")
    endif()
    execute_process(COMMAND sudo make install
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${OPENMPI_DIR}
                    RESULT_VARIABLE ompi_result
                    OUTPUT_VARIABLE ompi_ver)
    if (${hypre_result}==1)
        message("${hypre_ver}")
    endif()
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
    set(cmake_arg_config ../ -DHYPRE_INSTALL_PREFIX=${HYPRE_DIR})
    set(cmake_arg_build --build .)
    set(cmake_arg_install --install .)
    set(cmake_sudo sudo ${CMAKE_COMMAND})
    set(hypre_BUILD_dir _deps/hypre-src/src/cmbuild)
    execute_process(COMMAND ${CMAKE_COMMAND} ${cmake_arg_config}
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${hypre_BUILD_dir}
                    RESULT_VARIABLE hypre_result
                    OUTPUT_VARIABLE hypre_ver)
    if (${hypre_result}==1)
        message("${hypre_ver}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} ${cmake_arg_build}
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${hypre_BUILD_dir}
                    RESULT_VARIABLE hypre_result
                    OUTPUT_VARIABLE hypre_ver)
    if (${hypre_result}==1)
        message("${hypre_ver}")
    endif()
    execute_process(COMMAND ${cmake_sudo} ${cmake_arg_install}
                    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${hypre_BUILD_dir}
                    RESULT_VARIABLE hypre_result
                    OUTPUT_VARIABLE hypre_ver)
    if (${hypre_result}==1)
        message(FATAL_ERROR "${hypre_ver}")
    else()
        message("Hypre has been installed.")
        file(REMOVE_RECURSE build/_deps/hypre-src)
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
