#! /bin/bash

# run the setup script in the origin directory if it exists.
origin="/home/cclark/Code/sync/projects/HTMLSlideShowUtils"

if [ -x $origin/setup.sh ]
then
  echo "Running setup.sh script in $origin."
  $origin/setup.sh
  exit 1
fi

src=$(dirname $0)

echo -n "Copying utilities..."
cp $src/HTMLSlideShowUtils ./ -r
[ -d demo ] && rm -rf demo
cp $src/demo ./ -r
echo "DO NOT PUT ANYTHING IN THIS DIRECTORY YOU WANT TO KEEP" > demo/WARNING.txt
echo "It will be erased if the setup script is ever run again" >> demo/WARNING.txt
cp $src/setup.sh ./setup.sh
sed -i "/origin=/ s|/home/cclark/Code/sync/projects/HTMLSlideShowUtils|$src|" ./setup.sh
[ ! -e README.md ]   && cp HTMLSlideShowUtils/README.md ./
[ ! -e Makefile ]    && cp HTMLSlideShowUtils/Makefile ./
[ ! -e slides.md ]   && cp demo/slides.md ./
[ ! -e config.yaml ] && cp demo/config.yaml ./
echo "done"

echo -n "Checking dependencies..."

if ! hash pandoc > /dev/null 2>&1
then
  echo
  echo
  echo "'pandoc' is required, but is not installed."
  echo "Please install it before building your presentation."
  echo
fi


missing_modules=$(HTMLSlideShowUtils/scripts/check-imports.py HTMLSlideShowUtils/scripts/*.py | sort | uniq | grep -v 'macros')
if [ -n "${missing_modules}" ]
then
echo "These python modules are needed, but not installed"
echo "Please install them before building your presentation"
for module in ${missing_modules}
do
  echo "    ${module}"
done
fi
echo done
