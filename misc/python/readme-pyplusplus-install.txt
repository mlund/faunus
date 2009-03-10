website
-------

http://www.language-binding.net


gccxml
------

install
- ubuntu package - gccxml
- other distros - ?
- mac os - ?

attention, gccxml bug in ubuntu 8.10 intrepid ibex:
https://bugs.launchpad.net/ubuntu/+source/gccxml/+bug/293807
for now, fixed by the dummy mentioned in launchpad


pygccxml & py++
---------------

download both here:
http://sourceforge.net/project/showfiles.php?group_id=118209

unzip /home/andy/download/pygccxml-1.0.0.zip
cd pygccxml-1.0.0
sudo python setup.py install --prefix=/opt/pyplusplus-1.0.0
cd ..

unzip /home/andy/download/pyplusplus-1.0.0.zip
cd Py++-1.0.0
sudo python setup.py install --prefix=/opt/pyplusplus-1.0.0
cd ..

/opt/pyplusplus-1.0.0/env.sh:

case "$BASH_SOURCE" in
    /*)
    ENV_BASE_FILE="$BASH_SOURCE"
    ;;
    *)
    ENV_BASE_FILE="`pwd`/$BASH_SOURCE"
esac
ENV_BASE_DIR=`dirname $ENV_BASE_FILE`

export PATH=$ENV_BASE_DIR/bin:$PATH
export PYTHONPATH=$ENV_BASE_DIR/lib/python2.5/site-packages:$PYTHONPATH

unset BASE_FILE
unset ENV_BASE_DIR
