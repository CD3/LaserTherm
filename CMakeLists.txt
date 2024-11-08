cmake_minimum_required(VERSION 3.20)

project(LaserTherm)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/cmake/Modules")

include(externals/Plugins.cmake)

option(BUILD_UNIT_TESTS "Build unit tests for the library" ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${PROJECT_SOURCE_DIR}/cmake/Modules")

include(function-git_version)
git_version(LaserTherm)
message(STATUS "LaserTherm version: ${LaserTherm_VERSION}")
project(LaserTherm VERSION ${LaserTherm_VERSION})
# need to run configure_file twice
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/src/${PROJECT_NAME}/version.h.in"
               "${CMAKE_CURRENT_BINARY_DIR}/include/${PROJECT_NAME}/version.h")
configure_file("${CMAKE_CURRENT_BINARY_DIR}/include/${PROJECT_NAME}/version.h"
               "${CMAKE_CURRENT_BINARY_DIR}/include/${PROJECT_NAME}/version.h")

# dependencies
find_package(Boost REQUIRED)
find_package(Eigen3 3.3 REQUIRED)
find_package(libField REQUIRED)
find_package(UnitConvert REQUIRED)
find_package(OpenMP)


add_library(LaserTherm)
add_library(LaserTherm::LaserTherm ALIAS LaserTherm)
target_sources(
  LaserTherm
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/library.cpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Configuration/ptree_utils.cpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Configuration/Manager.cpp>
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/MaterialStructure.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/HeatSources/_1D/Cartesian/BeersLaw.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/HeatSources/_1D/FiniteDifferenceHeatSource.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Configuration/Builders.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Configuration/io.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Configuration/Manager.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Configuration/Units.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Configuration/ptree_utils.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Simulations/SingleEmitterExposure.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/HeatSolvers/BoundaryConditions.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/HeatSolvers/_1D/Cartesian/CrankNicholson.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/HeatSolvers/_1D/FiniteDifferenceHeatSolver.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/HeatSolvers/_1D/Cartesian/Explicit.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Materials/Types.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Materials/Basic.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Materials/Mixins.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Structures/_1D/AnyStructure.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Structures/_1D/Cartesian/Slab.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Structures/_1D/Infinite.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Utils/Concepts.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Utils/TypeTraits.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Utils/FiniteDifference.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Utils/TriDiagonalSolver.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Utils/Mixins.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Waveforms/ContinuousWave.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Waveforms/SinglePulse.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Waveforms/RegularPulseTrain.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Waveforms/Mixins.hpp>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/LaserTherm/Emitters/Basic.hpp>
)

target_include_directories(
  LaserTherm
  INTERFACE
    $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${${PROJECT_NAME}_BINARY_DIR}/include>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(
  LaserTherm
  INTERFACE Boost::boost Eigen3::Eigen libField::Field
    $<$<BOOL:OPENMP_CXX_FOUND>:OpenMP::OpenMP_CXX>
  PUBLIC UnitConvert::UnitConvert)

target_compile_definitions(LaserTherm
                           PUBLIC UNITCONVERT_NO_BACKWARD_COMPATIBLE_NAMESPACE)

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR}
   AND BUILD_UNIT_TESTS)
  enable_testing()
  add_subdirectory(testing)
endif()

# Install and Packaging
install(
  TARGETS LaserTherm
  EXPORT ${PROJECT_NAME}Targets
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  INCLUDES
  DESTINATION include/${PROJECT_NAME})
install(
  DIRECTORY src/
  DESTINATION include/${PROJECT_NAME}
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h")
install(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/
  DESTINATION include/${PROJECT_NAME}
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h")
install(
  EXPORT ${PROJECT_NAME}Targets
  FILE ${PROJECT_NAME}Targets.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION cmake)
file(
  WRITE ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  "include(CMakeFindDependencyMacro)
find_dependency(Boost)
find_dependency( Eigen3 )
include(\${CMAKE_CURRENT_LIST_DIR}/${PROJECT_NAME}Targets.cmake)
")
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${${PROJECT_NAME}_VERSION}
  COMPATIBILITY SameMajorVersion)
install(
  FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  DESTINATION cmake)

include(InstallRequiredSystemLibraries)

set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "")
set(CPACK_PACKAGE_VENDOR "")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.md")
if(${PROJECT_NAME}_VERSION_MAJOR)
  set(CPACK_PACKAGE_VERSION_MAJOR ${${PROJECT_NAME}_VERSION_MAJOR})
endif()
if(${PROJECT_NAME}_VERSION_MINOR)
  set(CPACK_PACKAGE_VERSION_MINOR ${${PROJECT_NAME}_VERSION_MINOR})
endif()
if(${PROJECT_NAME}_VERSION_PATCH)
  set(CPACK_PACKAGE_VERSION_PATCH ${${PROJECT_NAME}_VERSION_PATCH})
endif()
include(CPack)
