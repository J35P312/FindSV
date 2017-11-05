python setup.py
source FindSV_env.sh
cd SVDB
pip install -e .
cd ..
cd TIDDIT
chmod +x INSTALL.sh
./INSTALL.sh
cd ..
