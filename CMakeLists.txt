cmake_minimum_required(VERSION 3.17)

project(prony LANGUAGES C)
enable_language(C)

set(CMAKE_PREFIX_PATH D:/Apps/IntelSWTools/compilers_and_libraries/windows/mkl)
set(BLA_VENDOR Intel10_64lp_seq)
find_package(LAPACK)

add_executable(prony main.c)

# Intel MKL LAPACK and BLAS
target_include_directories(prony SYSTEM
    INTERFACE D:/Apps/IntelSWTools/compilers_and_libraries/windows/mkl/include)
include_directories(D:/Apps/IntelSWTools/compilers_and_libraries/windows/mkl/include)
target_link_directories(prony
    PUBLIC D:/Apps/IntelSWTools/compilers_and_libraries/windows/mkl/lib/intel64
    PUBLIC D:/Apps/IntelSWTools/compilers_and_libraries/windows/compiler/lib/intel64)
#link_libraries(mkl_intel_lp64_dll.lib mkl_intel_thread_dll.lib mkl_core_dll.lib libiomp5md.lib)
target_link_libraries(prony mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib)
link_libraries(mkl_intel_lp64_dll.lib mkl_sequential_dll.lib mkl_core_dll.lib)
#target_sources(prony
#    INTERFACE mkl_lapacke.h)
#target_sources(prony
#    PRIVATE lapacke.h
#    PRIVATE lapack.h
#    PRIVATE lapacke_mangling.h)

target_sources(prony
    PRIVATE main.c
    PUBLIC config.h
    PUBLIC types.h
    PRIVATE aux_print.c
    PUBLIC aux_print.h
    PRIVATE tclague.c
    PUBLIC tclague.h)
#target_link_libraries(prony
#    LINK_PUBLIC "${LAPACK_LIBRARIES}"
#    "${BLAS_mkl_intel_lp64_dll_LIBRARY}"
#    "${BLAS_mkl_core_dll_LIBRARY}"
#    "${BLAS_mkl_sequential_dll_LIBRARY}"
#)
