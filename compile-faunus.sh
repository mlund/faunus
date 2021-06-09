# INSTALL REQUIRED PYTHON PACKAGES
#conda create --name faunus 
#conda activate faunus
#conda install ruamel_yaml jsonschema jinja2 

source /opt/intel/oneapi/setvars.sh > /dev/null 2>&1

git checkout intel-clang
#git pull

export CXX=icpx
export CC=clang
#TBB_DIR="/glob/development-tools/versions/oneapi/2021.2/inteloneapi/tbb/2021.2.0/lib/cmake/tbb"
#cmake -DENABLE_OPENMP=off -DENABLE_TBB=on -DTBB_DIR=$TBB_DIR -DENABLE_CACHE=off -DPSTL_USE_PARALLEL_POLICIES=0
cmake -DENABLE_OPENMP=off -DENABLE_TBB=off -DENABLE_CACHE=off
make faunus -j

ctest -R unittests
