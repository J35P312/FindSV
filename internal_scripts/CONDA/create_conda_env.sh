cd $1
conda config --add channels r
conda config --add channels bioconda

conda env create --force --name FindSV_env -f $1/GENMOD.yml
