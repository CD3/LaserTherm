#! /usr/local/bin/python3
import sys
import subprocess as sp
import os
import numpy as np

# rebuild everything
	# create a new build2 directory
# make
# get output
# check for 'Error'/'H5Cpp.h'/etc
# delete new build2 directory
# exit based on if contained error or not

os.system("mkdir build2")
os.chdir("build2")
os.system("conan install .. --build missing")
os.system("cmake ..")
#ErrorMessage = os.popen("cmake --build .").read()
try:
    ErrorMessage = sp.run("cmake --build .", shell=True, check=True, text=True)
    os.chdir("..")
    os.system("rm -r build2")
    sys.exit(0)
except(sp.CalledProcessError):
    os.chdir("..")
    os.system("rm -r build2")
    sys.exit(1)
