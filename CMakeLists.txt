# BAM CMake file
# Humphrey Hu, 2015

cmake_minimum_required(VERSION 2.8)
project(bam)

if (NOT DEFINED CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif ()
message(STATUS "Build configuration: ${CMAKE_BUILD_TYPE}")

# Set version
set(bam_VERSION_MAJOR 0)
set(bam_VERSION_MINOR 2)
set(bam_VERSION ${bam_VERSION_MAJOR}.${bam_VERSION_MINOR})

# Set install paths
set (LIBRARY_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/lib" CACHE PATH
  "Target directory for all libraries.")
set (EXECUTABLE_OUTPUT_PATH "${PROJECT_SOURCE_DIR}/bin" CACHE PATH
  "Target directory for all executables.")

#set the debug flag and enable all warnings
add_definitions("-g -Wall -std=c++0x")

set(INCLUDE_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/include/bam")
set(LIBRARY_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(BINARY_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/bin")

include_directories( ${PROJECT_SOURCE_DIR}/include )

# TODO Set bit-sizes
add_library(bam SHARED src/bam.cc src/bam32Lookups.h)

set_target_properties(bam PROPERTIES
    PUBLIC_HEADER "include/bam.h")

add_executable(BAMTest tests/BAMTest.cc)
target_link_libraries(BAMTest bam)

install(TARGETS bam
    EXPORT bamTargets
    LIBRARY DESTINATION "${LIBRARY_INSTALL_DIR}" COMPONENT shlib
    PUBLIC_HEADER DESTINATION "${INCLUDE_INSTALL_DIR}"
    COMPONENT dev)

export(TARGETS bam
    FILE "${PROJECT_BINARY_DIR}/bamTargets.cmake")
    
configure_file(bamConfig.cmake.in
    "${PROJECT_BINARY_DIR}/bamConfig.cmake" @ONLY)

configure_file(bamConfigVersion.cmake.in
    "${PROJECT_BINARY_DIR}/bamConfigVersion.cmake" @ONLY)
    
install(FILES
    "${PROJECT_BINARY_DIR}/bamConfig.cmake"
    "${PROJECT_BINARY_DIR}/bamConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_PREFIX}" COMPONENT dev)
    
install(EXPORT bamTargets
    DESTINATION "${CMAKE_INSTALL_PREFIX}" COMPONENT dev)
