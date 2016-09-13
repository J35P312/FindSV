cd $1
conda config --add channels r
conda config --add channels bioconda

conda env create --force --name FindSV_env -f GENMOD.yml
conda install -y --name FindSV_env abyss
