cmake_minimum_required(VERSION 2.6)

project(minicsp)

option(BUILD_SHARED_LIBS OFF "True for building shared object")
set(CMAKE_CXX_FLAGS "-std=c++1z")
include_directories("/usr/local/include/")
include_directories(/usr/local/opt/libxml2/include/libxml2/)

find_package(ZLIB REQUIRED)

set(main "minicsp/xcsp3/main-xcsp3.cpp")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/minicsp/xcsp3")
# Basic Library
file(GLOB  SRC_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} minicsp/mtl/*.cpp minicsp/core/*.cpp minicsp/xcsp3/*.cpp)
#file(GLOB_RECURSE SRC_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
#list(REMOVE_ITEM lib_srcs ${main})

list(REMOVE_ITEM SRC_FILES minicsp/core/main.cpp)

include(ExternalProject)
ExternalProject_Add(XCSP3-CPP-Parser
        PREFIX parser
        GIT_REPOSITORY https://github.com/xcsp3team/XCSP3-CPP-Parser.git
        CONFIGURE_COMMAND ""
        BINARY_DIR parser/src/XCSP3-CPP-Parser/samples
        BUILD_COMMAND make lib
        INSTALL_COMMAND ""
        )

add_library(parserxcsp3core STATIC IMPORTED)
set_target_properties(parserxcsp3core PROPERTIES IMPORTED_LOCATION
        parser/src/XCSP3-CPP-Parser/lib/libparserxcsp3core.a)



if(BUILD_SHARED_LIBS)
    set(lib_type "SHARED")
else()
    set(lib_type "STATIC")
endif()


include_directories(
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${ZLIB_INCLUDE_DIR}
        #"/Users/audemard/csp/XCSP3-CPP-Parser/include/"
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake-build-debug/parser/src/XCSP3-CPP-Parser/include/
)

link_libraries(${ZLIB_LIBRARY})
link_libraries(parserxcsp3core)
link_libraries(xml2)

add_executable(minicsp ${main} ${SRC_FILES})
add_dependencies(minicsp XCSP3-CPP-Parser)
#target_link_libraries(minicsp)

