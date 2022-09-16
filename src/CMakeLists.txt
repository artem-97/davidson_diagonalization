get_filename_component(LIB_PATH "." ABSOLUTE)
get_filename_component(LIB_INCLUDE_PATH ".." ABSOLUTE)

file(GLOB_RECURSE LIB_CXX_SOURCES ${LIB_PATH}/*.cpp)
file(GLOB_RECURSE LIB_HEADERS ${LIB_PATH}/*.hpp)

add_library(davdiag ${LIB_CXX_SOURCES} ${LIB_HEADERS})

target_include_directories(davdiag PUBLIC ${LIB_INCLUDE_PATH})
target_link_libraries(davdiag PUBLIC pthread)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# find_package(BLAS)
# find_package(LAPACK)
target_link_libraries(davdiag PUBLIC Eigen3::Eigen)

# target_link_libraries(basis lapack)
if(CLANG_FORMAT_TOOL)
    add_clang_format_target(
        asmf_clang_format
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${LIB_HEADERS} ${LIB_CXX_SOURCES})
endif()

if(CLANG_TIDY_TOOL)
    add_clang_tidy_target(
        asmf_clang_tidy
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${LIB_INCLUDE_PATH}
        ${LIB_HEADERS} ${LIB_CXX_SOURCES})
endif()
