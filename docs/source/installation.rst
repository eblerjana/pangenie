==============
Installation
==============

We provide binaries for the latest version of PanGenie under `Releases <https://github.com/eblerjana/pangenie/releases/tag/v4.2.1>`_. Alternatively, PanGenie can be installed with conda::

  conda install bioconda::pangenie

However, if you prefer to build PanGenie yourself, follow the instructions below.

------------------------------------------------------
Building from source using Singularity (recommended)
------------------------------------------------------

Use the Singularity definition file located `here <https://github.com/eblerjana/pangenie/tree/master/container>`_ to build an (Ubuntu-based) container as follows (requires root privileges)::

 [sudo] singularity build pangenie.sif pangenie.def

In all usage examples below, call the ``PanGenie`` executables as follows::

 singularity exec pangenie.sif PanGenie-index <PARAMETERS>
 singularity exec pangenie.sif PanGenie <PARAMETERS>


For example, to show ``PanGenie``'s command line help, use the following command::

 singularity exec pangenie.sif PanGenie --help


You can check which versions of ``PanGenie`` (git hash) and of the ``jellyfish`` library have been installed in the container by running the following commands::

 singularity exec pangenie.sif cat /metadata/jellyfish.lib.version

should produce a line like this (so, here, v2.3.0)::

 $ libjellyfish-2.0-2:amd64 2.3.0-4build1 libjellyfish-2.0-dev:amd64 2.3.0-4build1

 singularity exec pangenie.sif cat /metadata/pangenie.git.version

should produce a line like this::

 $ 5a1f9c5


----------------------------------
Building from source using Conda
----------------------------------

::

 git clone https://github.com/eblerjana/pangenie.git
 cd pangenie
 conda env create -f environment.yml
 conda activate pangenie
 mkdir build; cd build; cmake .. ; make


-----------------------------------------------------------
Building from source
-----------------------------------------------------------

This requires jellyfish and cereal to be installed.

::

 git clone https://github.com/eblerjana/pangenie.git
 cd pangenie
 mkdir build; cd build; cmake .. ; make
