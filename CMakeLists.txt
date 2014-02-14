cmake_minimum_required(VERSION 2.8)
project(libbam)

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif ()

set(libbam_VERSION_MAJOR 0)
set(libbam_VERSION_MINOR 1)
set(libbam_VERSION ${libbam_VERSION_MAJOR}.${libbam_VERSION_MINOR})

#set the debug flag and enable all warnings
set(CMAKE_CXX_FLAGS "-g -Wall -std=c++0x")
set(INSTALL_CMAKE_DIR "${CMAKE_INSTALL_PREFIX}/lib/cmake/libbam" CACHE PATH
    "Installation directory for CMake files")

include_directories(
    ${libbam_SOURCE_DIR}/include
)

add_library(bam SHARED src/bam.cc src/bam32Lookups.h)

set_target_properties(bam PROPERTIES
    PUBLIC_HEADER "include/bam.h")

add_executable(BAMTest tests/BAMTest.cc)
target_link_libraries(BAMTest bam)

install(TARGETS bam
    EXPORT libbamTargets
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include/libbam"
)

export(TARGETS bam
    FILE "${PROJECT_BINARY_DIR}/libbamTargets.cmake")
export(PACKAGE libbam)

set(CONF_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(CONF_LIBRARY "${CMAKE_INSTALL_PREFIX}/lib/libbam.so")

install(EXPORT libbamTargets
    DESTINATION "${INSTALL_CMAKE_DIR}")

# #Lets build the documentation
# #######################################################################################
# #-- Add an Option to toggle the generation of the API documentation
# #option(BUILD_DOCUMENTATION "Use Doxygen to create the HTML based API documentation" ON)
# #if(BUILD_DOCUMENTATION)
#   FIND_PACKAGE(Doxygen)
#   if (NOT DOXYGEN_FOUND)
#     message(FATAL_ERROR 
#       "Doxygen is needed to build the documentation. Please install it correctly")
#   endif()
#   #-- Configure the Template Doxyfile for our specific project
#   configure_file(Doxyfile ${localizer_SOURCE_DIR}/Doxyfile  @ONLY IMMEDIATE)
#   #-- Add a custom target to run Doxygen when ever the project is built
#   add_custom_target (docs #ALL 
#     COMMAND ${DOXYGEN_EXECUTABLE} ${localizer_SOURCE_DIR}/Doxyfile
#     SOURCES ${localizer_SOURCE_DIR}/Doxyfile)
# 
#   # IF you do NOT want the documentation to be generated EVERY time you build the project
#   # then leave out the 'ALL' keyword from the above command.
# #endif()
# ##############################################################################################
