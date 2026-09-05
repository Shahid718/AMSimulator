# Only for Linux, nvfortran
# Modified : 09 September 2026

cmake_minimum_required(VERSION 3.12)
project(AMSimulator VERSION 1.0 LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

set(PARALLELIZATION "OPENMP" CACHE STRING "Choose parallelization method: OPENMP or DO_CONCURRENT")
set_property(CACHE PARALLELIZATION PROPERTY STRINGS OPENMP DO_CONCURRENT)

if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(STATUS "Detected Linux operating system")
    set(LINUX TRUE)
else()
    message(FATAL_ERROR "This project is designed for Linux OS only. Detected: ${CMAKE_SYSTEM_NAME}")
endif()

# --- Compiler configuration ---
# Currently only tested against NVIDIA HPC SDK (nvfortran). Other compilers
# may work but are untested and unsupported until verified.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")
    message(STATUS "Using NVIDIA Fortran compiler: ${CMAKE_Fortran_COMPILER}")
    set(COMPILER_NVFORTAN TRUE)

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-g -O0 -Mbounds")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -fast")

    set(OpenMP_Fortran_FLAGS "-mp")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,stack-size=500000000")
    set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

else()
    message(FATAL_ERROR
        "Unsupported or untested compiler: ${CMAKE_Fortran_COMPILER_ID}.\n"
        "This project is currently only tested with the NVIDIA HPC SDK (nvfortran).\n"
        "Support for other compilers (e.g. Intel) may be added in a future release "
        "once verified. Contributions/testing reports are welcome.")
endif()

# --- Parallelization backend ---
if(PARALLELIZATION STREQUAL "OPENMP")
    message(STATUS "Using OpenMP parallelization")

    find_package(OpenMP REQUIRED)
    if(OpenMP_Fortran_FOUND)
        message(STATUS "OpenMP found")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
        add_compile_definitions(OPENMP)
    else()
        message(FATAL_ERROR "OpenMP not found. Please install OpenMP or use DO_CONCURRENT")
    endif()

elseif(PARALLELIZATION STREQUAL "DO_CONCURRENT")
    message(STATUS "Using DO CONCURRENT parallelization")
    add_compile_definitions(DC)

    if(COMPILER_NVFORTAN)
        # Without -stdpar, nvfortran compiles DO CONCURRENT as a single-threaded
        # (SIMD-vectorized) loop. This flag is required for actual CPU multithreading.
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -stdpar=multicore")
    endif()

else()
    message(FATAL_ERROR "Invalid PARALLELIZATION option: ${PARALLELIZATION}. Choose OPENMP or DO_CONCURRENT")
endif()

# --- Output directories ---
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/output)

set(SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/PrecisionRangeModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/InputManagerModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ErrorModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/TimeStampModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/InitialMicrostructureModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/AnisotropySolverModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/DrivingForceModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/StencilModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/UpdateFieldModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/NoiseModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/OutputManagerModule.f90
    ${CMAKE_CURRENT_SOURCE_DIR}/src/ArrayManagerModule.f90
)

set(MAIN_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/app/main.f90)

add_executable(main ${SOURCES} ${MAIN_SOURCE})

set_target_properties(main PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
)

target_include_directories(main PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_SOURCE_DIR}/app
)

add_custom_command(TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_BINARY_DIR}/main
        ${CMAKE_BINARY_DIR}/output/
    COMMENT "Copying executable to output directory"
)

add_custom_command(TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo "============================================"
    COMMAND ${CMAKE_COMMAND} -E echo "Build completed successfully!"
    COMMAND ${CMAKE_COMMAND} -E echo "Executable: ${CMAKE_BINARY_DIR}/output/main"
    COMMAND ${CMAKE_COMMAND} -E echo "Compiler: ${CMAKE_Fortran_COMPILER_ID}"
    COMMAND ${CMAKE_COMMAND} -E echo "Parallelization: ${PARALLELIZATION}"
    COMMAND ${CMAKE_COMMAND} -E echo "Build type: ${CMAKE_BUILD_TYPE}"
    COMMAND ${CMAKE_COMMAND} -E echo "============================================"
    COMMENT "Build information"
)

message(STATUS "Fortran flags: ${CMAKE_Fortran_FLAGS}")
message(STATUS "Fortran flags Debug: ${CMAKE_Fortran_FLAGS_DEBUG}")
message(STATUS "Fortran flags Release: ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "Linker flags: ${CMAKE_EXE_LINKER_FLAGS}")
message(STATUS "Parallelization: ${PARALLELIZATION}")
if(PARALLELIZATION STREQUAL "OPENMP")
    message(STATUS "OpenMP flags: ${OpenMP_Fortran_FLAGS}")
endif()