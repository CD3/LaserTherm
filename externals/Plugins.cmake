if(EXISTS ${CMAKE_CURRENT_LIST_DIR}/Format.cmake)
  message(STATUS "Adding Format.cmake")
  set(FORMAT_CHECK_CMAKE ON)
  set(CMAKE_FORMAT_EXCLUDE
      ".*/Format.cmake/.*|.*/Format.cmake"
      CACHE STRING "CMake formatting file exclusion pattern.")
  add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/Format.cmake)
else()
  message(STATUS "Format.cmake not found..")
endif()
