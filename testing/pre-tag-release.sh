#! /bin/bash

set -e
root=$(git rev-parse --show-toplevel)

cd $root

bindir="$PWD/build-and-test"
function cleanup()
{
  rm -r $bindir
}
trap cleanup EXIT

function copy_bindir()
{
  cp -r $bindir $bindir.error
}
trap copy_bindir ERR

echo "Checking that project builds."
mkdir $bindir
cd $bindir
conan install ${root} -g cmake_paths
cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_paths.cmake -DCMAKE_INSTALL_PREFIX=$bindir/install
cmake --build .

echo "Checking that unit tests pass."
cmake --build . --target test

echo "Checking that project can be installed."
cmake --build . --target install




echo "WARNING: test to check that the installed project can be detected and used by CMake have not been implemented yet."
mkdir app
cd app

cat << EOF > main.cpp
#include <iostream>

int main()
{

  return 0;
}
EOF

cat << EOF > CMakeLists.txt
cmake_minimum_required(VERSION 3.1)
add_executable( main main.cpp )
# find_package( TemplateProject PREQUIRED )
# target_link_libraries(main TemplateProject::TemplateProject )
set_target_properties(main PROPERTIES CXX_STANDARD 11)
EOF

mkdir build2
cd build2
conan install ${root}
cmake .. -DCMAKE_INSTALL_PREFIX=$bindir/install
cmake --build .
./main

cd ..







echo "PASSED"

exit 0
