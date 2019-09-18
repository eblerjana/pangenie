# PGG-Typer

## Requirements
* gcc 4.7+
* cmake
* jellyfish

## Installation
`` git clone https://jana_ebler@bitbucket.org/jana_ebler/pgg-typer.git``  
`` cd pgg-typer``  
``mkdir build; cd build; cmake .. ; make``

## Installing into a conda environment
`` git clone https://jana_ebler@bitbucket.org/jana_ebler/pgg-typer.git``  
`` cd pgg-typer``  
`` conda env create -f environment.yml``  
`` conda activate pggtyper``  
``export PKG_CONFIG_PATH="/MMCI/TM/scratch/jebler/miniconda3/envs/pggtyper/lib/pkgconfig"``  
``mkdir build; cd build; cmake .. ; make``

