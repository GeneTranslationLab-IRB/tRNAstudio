#!/bin/bash

# Check conda installation
conda="which conda"
condaInstalled=$(eval "$conda")

# Install conda 
if [[ "$condaInstalled" == "" ]]
then
    if [[ "$OSTYPE" == linux-gnu* ]]
    then
    #Download anaconda file for Linux
    wget https://repo.anaconda.com/archive/Anaconda3-2020.07-Linux-x86_64.sh
    #Checking the integrity of the file
    sha256sum Anaconda3-2020.07-Linux-x86_64.sh
    #Running the .sh script
    bash Anaconda3-2020.07-Linux-x86_64.sh
    #Compiling from source
    echo 'PATH=~/anaconda3/bin:$PATH' >> ~/.bashrc
    export PATH=~/anaconda3/bin:$PATH
    source ~/.bashrc
    rm Anaconda3-2020.07-Linux-x86_64.sh
    fi


    if [[ "$OSTYPE" = darwin* ]]
    then
    #Download Homebrew 
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
    #Download wget
    brew install wget
    #Download anaconda file for MAC
    wget https://repo.anaconda.com/archive/Anaconda3-2020.07-MacOSX-x86_64.sh
    #Checking the integrity of the file
    sha256sum Anaconda3-2020.07-MacOSX-x86_64.sh
    #Running the .sh script
    bash Anaconda3-2020.07-MacOSX-x86_64.sh
    #Compiling from source
    echo 'PATH=~/anaconda3/bin:$PATH' >> .bashrc
    export PATH=~/anaconda3/bin:$PATH
    source ~/.bashrc
    rm Anaconda3-2020.07-MacOSX-x86_64.sh
    fi

fi

# Create conda enviroment (tRNA-PipelineEnv)
conda env create --name tRNAstudioEnv --file=environment.yml





