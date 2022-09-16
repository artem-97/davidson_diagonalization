cmake_minimum_required(VERSION 3.14)

project(davdiag)

include(cmake/ClangFormat.cmake)
include(cmake/ClangTidy.cmake)
include(cmake/CompileOptions.cmake)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "bin")

add_subdirectory(src)

if(DAVDIAG_TESTS OR DAVDIAG_DEVELOPER)
    add_subdirectory(tests)
endif()
