cmake_minimum_required(VERSION 3.20)

find_package(Catch2 REQUIRED)

# Find all sources.
file(
  GLOB_RECURSE unitTest_SOURCES
  RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  "./CatchTests/*.cpp")
message(STATUS "Detected Catch-based Unit Test Sources:")
foreach(testSrc ${unitTest_SOURCES})
  message(STATUS "  ${testSrc}")
endforeach()

set(testName ${PROJECT_NAME}_CatchTests)
add_executable(${testName} ${unitTest_SOURCES})
target_compile_definitions(${testName} PRIVATE -DTESTING)
target_link_libraries(${testName} Catch2::Catch2WithMain LaserTherm)
target_include_directories(${testName}
                           PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
set_target_properties(${testName} PROPERTIES CXX_STANDARD 17)

add_test(
  NAME ${testName}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${testName})
