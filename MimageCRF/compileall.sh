# Use one of these values for OS_TYPE in the make commands below :
# linux32 or linux64 or macosx 

cd zlib-1.2.3/
make clean
make
cd ..

cd imageLib/
make allclean
make
cd ..

cd MRF2.0float/
make allclean
make OS_TYPE=linux64
cd ..

cd libconfig-1.4.6/
make clean
./configure
make
rm lib/.libs/libconfig++.so*
rm lib/.libs/libconfig*.dylib
cd ..

cd opticalflow/
make clean
make
cd ..

make clean
make OS_TYPE=linux64
