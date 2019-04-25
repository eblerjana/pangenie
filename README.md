# PGG-Typer

## Requirements
* gcc 4.7+
* cmake

## Installation
`` git clone https://jana_ebler@bitbucket.org/jana_ebler/pgg-typer.git``  
`` cd pgg-typer``  
``mkdir build; cd build; cmake .. ; make``

## Installing into a conda environment
`` git clone https://jana_ebler@bitbucket.org/jana_ebler/pgg-typer.git``  
`` cd pgg-typer``  
`` conda env create -f environments.yml``  
``export CFLAGS="-I/MMCI/TM/scratch/jebler/miniconda3/envs/pgg-typer/include"``  
``export LDFLAGS="-L/MMCI/TM/scratch/jebler/miniconda3/envs/pgg-typer/lib"``  
``export CPPFLAGS="-I/MMCI/TM/scratch/jebler/miniconda3/envs/pgg-typer/include"``  
``export CPATH="/MMCI/TM/scratch/jebler/miniconda3/envs/pgg-typer/include"``  
``mkdir build; cd build; cmake .. ; make``

