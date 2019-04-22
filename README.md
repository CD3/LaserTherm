# LaserTherm

## Description

`LaserTherm` is a C++ library, and collection of client applications, for running laser-tissue interaction simulations. At its
core, `LaserTherm` is a heat solver, with models to incorporate heat due to laser exposure. It was created with the following
design goals:

1. Simple
1. Modular
1. Flexible

## Building and Installing

Cmake is used to configure, build, and install:

```
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
$ cmake --build . --target install
```
