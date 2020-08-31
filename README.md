This archive contains the C++ implementation of our paper, which is integrated in a custom version of TTK (``ttk`` directory), as well as scripts to reproduce performance numbers (``benchmarks``, see the README.md there).

# Building and installing our implementation
The following instructions describe the installation process for a recent Ubuntu Linux operating system. These instructions may need to be adapted if you use a different version or a different operating system.

## Installing the dependencies
Several dependencies need to be installed. Under Ubuntu:

    $ sudo apt-get install cmake-qt-gui
    $ sudo apt-get install libvtk7-dev
    $ sudo apt-get install qt5-default qttools5-dev libqt5x11extras5-dev

## Downloading ParaView
Our implementation is based on TTK. It requires to build ParaView **v5.7.0** from source. Note that the code will not build for any other version. Under ubuntu:

    $ wget "https://www.paraview.org/paraview-downloads/download.php?submit=Download&version=v5.7&type=source&os=Sources&downloadFile=ParaView-v5.7.0.tar.gz" -O ParaView-v5.7.0.tar.gz
    $ tar xvzf ParaView-v5.7.0.tar.gz
    
## Patching ParaView
Now, go to our source directory (``ttk``), specifically under the subdirectory ``ttk/paraview/patch``. From there:

    $ ./patch-paraview-5.7.0.sh <path to your decompressed ParaView source tree, cf above>

## Configuring, building and installing ParaView
In ParaView's source tree, create a ``build`` subdirectory and enter it. There, call cmake as follows:

    $ cmake .. \
      -DCMAKE_BUILD_TYPE=Release \
      -DPARAVIEW_ENABLE_PYTHON=ON \
      -DPARAVIEW_INSTALL_DEVELOPMENT_FILES=ON \
      -DVTK_PYTHON_VERSION=3
      
Then, to build ParaView, enter the following command, where ``N`` is the number of available cores on your system (this will take a **LONG** time):

    $ make -jN
    
To install your build of ParaView, enter the following command (under Ubuntu Linux):

    $ sudo make install
    
## Configuring, building and installing our custom version of TTK
Now, go to our source directory (``ttk``) and create a ``build`` subdirectory and enter it. There, call cmake as follows:

    $ cmake .. 
      
Then, to build TTK, enter the following command, where ``N`` is the number of available cores on your system (this will take a **LONG** time):

    $ make -jN
    
To install our custom version of TTK, enter the following command (under Ubuntu Linux):

    $ sudo make install
    
# Reproducing the performance numbers 
To reproduce the performance figures from our paper, please enter the ``benchmarks`` directory of our archive, and see the README.md there for further instructions.
