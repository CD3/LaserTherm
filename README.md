# LaserTherm

## Description

`LaserTherm` is a C++ library, and collection of client applications, for running laser-tissue interaction simulations. At its
core, `LaserTherm` is a collection of heat solvers, with models to incorporate heat due to laser exposure. It was created with the following
design goals:

1. Simple
1. Modular
1. Flexible

It is currently being developed as a series of undergraduate research projects in the [Department of Physics](https://www.fhsu.edu/physics/) at [Fort Hays State University](https://www.fhsu.edu). 

Authors:

- Dr. CD Clark III
- Daniel Huantes
- Emily Cranwell

## Building and Installing

If not already installed, install python3, pip3, git, cmake, and the gnu C++ compiler. You will need administrator
privileges, or ask your system administrator to install them for you. It is likely
that they are already installed.

```
$ sudo apt install python3 python3-pip git cmake gcc g++
```
Everything else can be installed to your user account and does not require admin
privilege.

#### Install Conan
Conan is a C/C++ package manager written in python. You can use it to download, build, and install
all of the `LaserTherm` dependencies. Install Conan with `pip3`.
```
pip3 install conan --user
```
The `--user` option will cause Conan will be installed in a subdirectory named
`.local/` under your users home directory. To run `conan` like a regular command
you need to add `${HOME}/.local/bin` to your `PATH` variable.
Add the following lines,
```
PATH="${PATH}:${HOME}/.local/bin"
export PATH
```
to the end of your `.bashrc` file (Note the **double** quotes and the **colon**). It should be located in your home directory
(it will be hidden since it starts with a '.'), and it is OK to create it if it
doesn't.  For example, if you use `vim`,
```
$ vim ~/.bashrc
# add the lines above and save
```
You will need to open a new terminal for the change to take effect. After
that, you should be able to run `conan`
```
$ conan
Consumer commands
  install    Installs the requirements specified in a recipe (conanfile.py or conanfile.txt).
  config     Manages Conan configuration.
  get        Gets a file or list a directory of a given reference or package.
  info       Gets information about the dependency graph of a recipe.
  search     Searches package recipes and binaries in the local cache or in a remote.
Creator commands
  new        Creates a new package recipe template with a 'conanfile.py' and optionally,
             'test_package' testing files.
  create     Builds a binary package for a recipe (conanfile.py).
  upload     Uploads a recipe and binary packages to a remote.
  export     Copies the recipe (conanfile.py & associated files) to your local cache.
  export-pkg Exports a recipe, then creates a package from local source and build folders.
  test       Tests a package consuming it from a conanfile.py with a test() method.
Package development commands
  source     Calls your local conanfile.py 'source()' method.
  build      Calls your local conanfile.py 'build()' method.
  package    Calls your local conanfile.py 'package()' method.
  editable   Manages editable packages (package that resides in the user workspace, but are
             consumed as if they were in the cache).
  workspace  Manages a workspace (a set of packages consumed from the user workspace that
             belongs to the same project).
Misc commands
  profile    Lists profiles in the '.conan/profiles' folder, or shows profile details.
  remote     Manages the remote list and the package recipes associated to a remote.
  user       Authenticates against a remote with user/pass, caching the auth token.
  imports    Calls your local conanfile.py or conanfile.txt 'imports' method.
  copy       Copies conan recipes and packages to another user/channel.
  remove     Removes packages or binaries matching pattern from local cache or remote.
  alias      Creates and exports an 'alias package recipe'.
  download   Downloads recipe and binaries to the local cache, without using settings.
  inspect    Displays conanfile attributes, like name, version and options. Works locally, in
             local cache and remote.
  help       Shows help for a specific command.
  graph      Generates and manipulates lock files.

Conan commands. Type "conan <command> -h" for help
```
Note: if don't add the `pip3` local bin directory to your path, you can still run Conan directly
```
$ ~/.local/bin/conan
```
Adding it to your `PATH` just makes it available as a command.

You also need to set the compiler settings in conan to avoid certain version conflicts with the following comman
```
$ conan profile update settings.compiler.libcxx=libstd11 default
```

### Install `LaserTherm` Dependencies

You need to add Dr. Clark's Conan package repository to download the the `LaserTherm
```
$ conan remote add cd3 https://api.bintray.com/conan/cd3/conan-devel
```

Now download this repository
```
$ git clone http://fermi.fhsu.edu:81/CD3/Laser
Username for 'http://fermi.fhsu.edu:81':
Password for 'http://fermi.fhsu.edu:81':
```
You will have to enter your fermi username and password.

This repository contains git "submodules", which must also be downloaded. You do this by first
initializing, and then updating, the submodules.
```
$ git submoule init
$ git submoule update
```


### Building `LaserTherm`

Now you can configure and build. `LaserTherm` uses CMake, so first create a `build/`
directory to compile the library in. Go into this directory and use conan to download, build,
and install the dependencies. Conan will create a script named `activate.sh` that contains all
of the information needed for CMake to use the packages installed by Conan (in fact, we use
conan to install CMake because the version currently available in the Ubuntu repository is too old).
To setup your environment to run CMake, you need to source the `activate.sh` script.
You can then run `cmake` to configure the build, and then again to compile the build.
To do all of this, run the following commands.
```
$ cd LaserTherm
$ mkdir build
$ cd build
$ conan install .. --build missing
$ source activate.sh  # for Linux or Mac
$ cmake ..
$ cmake --build .
```
If everything went well, you should see a directory named `testBin/` that contains a file named
`LaserTherm`. You can now run the unit tests.
```
(conanenv) $ ./testBin/LaserTherm_CatchTests
```
Congratulations! Have fun.

