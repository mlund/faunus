# INSTALL REQUIRED PYTHON PACKAGES
#conda create --name faunus 
#conda activate faunus
#conda install ruamel_yaml jsonschema jinja2 

git checkout intel-clang
git pull

export CXX=icpx
export CC=clang
TBB_DIR="/glob/development-tools/versions/oneapi/2021.2/inteloneapi/tbb/2021.2.0/lib/cmake/tbb"
cmake -DENABLE_OPENMP=off -DENABLE_TBB=on -DTBB_DIR=$TBB_DIR
make faunus -j

ctest -R unittests
