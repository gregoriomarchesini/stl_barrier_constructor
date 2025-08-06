cmake_minimum_required(VERSION 3.10)
project(HypercubeLib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/lib/mosek/11.0/tools/platform/linux64x86/h

)

file(GLOB_RECURSE HYPERCUBE_SOURCES "${PROJECT_SOURCE_DIR}/src/*.cpp")
add_library(hypercube ${HYPERCUBE_SOURCES})

function(add_example name)
    add_executable(${name} examples/${name}.cpp)
    target_link_directories(${name} PRIVATE
        ${PROJECT_SOURCE_DIR}/lib/mosek/11.0/tools/platform/linux64x86/bin
    )
    target_link_libraries(${name} PRIVATE
        hypercube
        mosek64
        mosekxx11_0
        tbb
        fusion64
    )
    set_target_properties(${name} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
        INSTALL_RPATH_USE_LINK_PATH TRUE
    )
endfunction()

message(STATUS "🔧 Running 'make install' in Fusion C++ directory...")
execute_process(
    COMMAND make install
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/lib/mosek/11.0/tools/platform/linux64x86/src/fusion_cxx
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error_output
)



# Then just call it
add_example(minimal_ex)