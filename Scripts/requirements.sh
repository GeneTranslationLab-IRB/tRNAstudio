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
    source ~/.bashrc
    rm Anaconda3-2020.07-Linux-x86_64.sh
    fi


    if [[ "$OSTYPE" = darwin* ]]
    then
    #Download Homebrew 
    /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)"
    #Download anaconda and miniconda 
    brew install --cask anaconda
    brew cask install miniconda
    #install wget
    fi

fi

# Create conda enviroment (tRNA-PipelineEnv)
conda env create --name tRNA-PipelineEnv --file=environment.yml





