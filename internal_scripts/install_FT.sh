cd $1
git clone https://github.com/J35P312/TIDDIT.git
cd TIDDIT
git checkout -b J35P312 2.0
mkdir build
cd build
cmake ..
make
