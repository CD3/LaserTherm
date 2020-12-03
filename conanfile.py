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
    requires = 'boost/1.72.0','eigen/3.3.7@cd3/devel', 'libField/0.8@cd3/devel', 'UnitConvert/0.11@cd3/devel'

    def build(self):

      cmake = CMake(self)

      # get location of the eigen cmake config file to pass to cmake
      defs = {}
      cmake.configure(defs=defs)

      cmake.build()
