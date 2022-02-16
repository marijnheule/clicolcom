# need ninja, meson, boost
cd ChromSAT
ln -s ../minicsp/minicsp
ln -s ../sparsehash/src/sparsehash
ln -s ../tclap/include/tclap
#ln -s /opt/homebrew/Cellar/boost/1.76.0/include/boost
cd ../sparsehash
./configure; make; make install
cd ../ChromSAT
meson --buildtype=debugoptimized build
meson compile -C build
