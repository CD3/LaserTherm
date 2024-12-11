list:
  just --list

install-deps:
  conan install . --build missing -s build_type=Debug
  conan install . --build missing -s build_type=Release

configure:
  cd build && cmake -DCMAKE_TOOLCHAIN_FILE="generators/conan_toolchain.cmake" ..

build-debug:
  cmake --build build --config Debug

build-release:
  cmake --build Release

run-tests:
  cd build && ./testing/Debug/LaserTherm_CatchTests
