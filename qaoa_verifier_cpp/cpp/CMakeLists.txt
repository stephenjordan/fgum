cmake_minimum_required(VERSION 3.16)
project(qaoa_verifier CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMP: optional (AppleClang doesn't ship it by default)
find_package(OpenMP QUIET)

# Float64 sources: aggressive optimization
set(F64_SOURCES
    primitives.cpp
    branch.cpp
    root.cpp
    contract.cpp
    alloc.cpp
)
add_library(qaoa_f64 OBJECT ${F64_SOURCES})
target_compile_options(qaoa_f64 PRIVATE -O3)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(qaoa_f64 PRIVATE -march=native)
endif()
if(OpenMP_CXX_FOUND)
    target_link_libraries(qaoa_f64 PRIVATE OpenMP::OpenMP_CXX)
endif()

# DD sources: strict IEEE 754 (no FMA for Dekker splitting correctness)
set(DD_SOURCES
    dd.cpp
    dd_fused.cpp
)
add_library(qaoa_dd OBJECT ${DD_SOURCES})
target_compile_options(qaoa_dd PRIVATE -O3 -ffp-contract=off)
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    target_compile_options(qaoa_dd PRIVATE -march=native)
    include(CheckCXXCompilerFlag)
    check_cxx_compiler_flag(-mno-fma HAS_MNO_FMA)
    if(HAS_MNO_FMA)
        target_compile_options(qaoa_dd PRIVATE -mno-fma)
    endif()
endif()
if(OpenMP_CXX_FOUND)
    target_link_libraries(qaoa_dd PRIVATE OpenMP::OpenMP_CXX)
endif()

# Verifier executable
add_executable(verify
    verify.cpp
    $<TARGET_OBJECTS:qaoa_f64>
    $<TARGET_OBJECTS:qaoa_dd>
)
if(OpenMP_CXX_FOUND)
    target_link_libraries(verify PRIVATE OpenMP::OpenMP_CXX)
endif()

# Shared library (optional, for Python ctypes if needed)
add_library(qaoa_contract_shared SHARED
    exports.cpp
    $<TARGET_OBJECTS:qaoa_f64>
    $<TARGET_OBJECTS:qaoa_dd>
)
set_target_properties(qaoa_contract_shared PROPERTIES
    OUTPUT_NAME qaoa_contract
    POSITION_INDEPENDENT_CODE ON
)
if(OpenMP_CXX_FOUND)
    target_link_libraries(qaoa_contract_shared PRIVATE OpenMP::OpenMP_CXX)
endif()

# PIC for object libraries (needed by shared lib)
set_target_properties(qaoa_f64 PROPERTIES POSITION_INDEPENDENT_CODE ON)
set_target_properties(qaoa_dd PROPERTIES POSITION_INDEPENDENT_CODE ON)
