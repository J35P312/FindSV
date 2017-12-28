if [ "$1" == "noinstall" ]
    then
        python setup.py 0
elif [ "$1" == "" ]
    then

        python setup.py 1 
        source FindSV_env.sh
        cd SVDB
        pip install -e .
        cd ..
        cd TIDDIT
        chmod +x INSTALL.sh
        ./INSTALL.sh
        cd ..

else
        echo "invalid option" 
        echo "./setup.sh                  --- install software and set the path to reference, scripts, etc"
        echo "./setup.sh noinstall        --- skip the installation of software (these needs to be installed manually)"

    fi



