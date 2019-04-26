from conans import ConanFile, CMake, tools
import os

class ConanBuild(ConanFile):
    name = "LaserTherm"
    version = "master"
    license = "MIT"
    author = "CD Clark III clifton.clark@gmail.com"
    url = "http:/fermi.fhsu.edu/CD3/LaserTherm"
    homepage = "http://fermi.fhsu.edu/CD3/LaserTherm"
    description = "A C++ library for doing laser-tissue simulation.."
    topics = ("C++", "Laser", "Physics")

    generators = "cmake", "virtualenv"
    requires = 'boost/1.69.0@conan/stable','eigen/3.3.7@conan/stable', 'libField/master@local/devel', 'UnitConvert/master@local/devel'
    build_requires = 'cmake_installer/3.13.0@conan/stable', 'ninja_installer/1.9.0@bincrafters/stable'

    def build(self):

      cmake = CMake(self)

      # get location of the eigen cmake config file to pass to cmake
      defs = {}
      defs["Eigen3_DIR"] = os.path.join( self.deps_cpp_info["eigen"].rootpath, "share", "eigen3", "cmake" )
      cmake.configure(defs=defs)

      cmake.build()
