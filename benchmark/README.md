This file contains instructions to reproduce the results from the paper. 
We assume that readers are using a Unix system with Bash. We also assume that they successfully installed our modified version of TTK containing the implementation of our paper (available as additional material) along with ParaView 5.7.0.

Before following the instructions below, please enter the `scripts' sub-directory.

# 0. Setup
To setup the following scripts, please edit the file `setup.sh' to adjust the necessary variables. In particular, you will need to adjust the following environment variables:
  - `PVPYTHON_DIR': directory where the binary `pvpython' has been installed
    (default value: /usr/local/bin/)
  - `PV_PLUGIN_PATH': directory where the TTK ParaView plugins have been
    installed (default value: /usr/local/bin/plugins/)

# 1. Data download
The file `url.txt' contains the URL of each data set used in this paper.
To download them all in the current directory, use the script `download.sh':
<code>
./download.sh
</code>

# 2. Convert data to VTK file format
Next, to convert the downloaded files to the VTK file format, you can use the script `convert.sh':
<code>
./convert.sh
</code>

# 3. Generate the random data set
Next, to generate the random data set, you can use the script `random.sh':
<code>
./random.sh
</code>

# 4. Pre-process the oceanVortices data set
The data set `oceanVortices.vtu' needs some pre-processing in order to assign data values on vertices instead of faces and to isolate continents. To run this pre-processing, run the script `preProcessOceanVortices.sh` as follows:
<code>
./preProcessOceanVortices.sh
</code>

# 5. Reproducing the first performance table
To reproduce the first performance table (Table 1, documenting the running times of the simplification), run the following script (this will take a LONG time):
<code>
./performance_simplification.sh
</code>
Next, to convert the raw output of that script into a LaTeX table, run the following script (this will produce the file `performance_simplification.tex'):
<code>
./performance_simplification_2latex.sh
</code>

# 6. Reproduce the random stress curves
To reproduce the curves documenting the random stress experiment (Figure 7, time as a function of the percentage of removed extrema), run the following script (this will take a LONG time):
<code>
./performance_randomSimplification.sh
</code>
Next, to convert the raw output of that script into a LaTeX curve, run the following script (requires `gnuplot', this will produce the file `plotRandomSimplification.tex'):
<code>
./plot_randomSimplification.gplot
</code>

# 7. Reproduce the second performance table
To reproduce the second performance table (Table 2, documenting the running times for the persistence-driven simplification), run the following script (this will take a LONG time):
<code>
./performance_persistenceSimplification.sh
</code>
Next, to convert the raw output of that script into a LaTeX table, run the following script (this will produce the file `performance_persistenceSimplification.tex'):
<code>
./performance_persistenceSimplification_2latex.sh
</code>
