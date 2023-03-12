LAIpath: Path Length Distribution Model for LAI estimation 

The path-length distribution method (PATH) was proposed
by Hu et al. (2014) to address not only the clumping effect
between canopies but also the crown-shape-induced nonrandomness
within canopies. It was developed based on the path length
distribution and a theoretical prototype of Beer¡¯s law.
Path length distribution theory can characterize and handle
nonrandomness within canopies, which may cause an underestimation
of LAI of up to 25%. Previous algorithms cannot deal
with such a problem.

Detailed theory is discribed in:
1. Hu, R. et al. (2014). Indirect measurement of leaf area index on the basis of path length distribution. REMOTE SENSING OF ENVIRONMENT, 155, 239-247. doi: http://dx.doi.org/10.1016/j.rse.2014.08.032

2. G. Yan, R. Hu, J. Luo, M. Weiss, H. Jiang, X. Mu, D. Xie, W. Zhang et al. (2019). Review of indirect optical measurements of leaf area index: Recent advances, challenges, and perspectives. Agricultural and Forest Meteorology, 265, 390-411. doi:http://dx.doi.org/10.1016/j.agrformet.2018.11.033

It is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License. Please cite the papers above if you use it for publishing.

The codes are based on GNU Scientific Library (GSL).
It is convenient to install and use GSL under Linux.
Visual Studio project file and necessary libraries are provided under Windows.

1. An example for using Path Length Distribution Model
example.cpp

2. Soucre code for Path Length Distribution Model
LAIPath.cpp 
LAIPath.h  

3. Visual Studio project file 
LAI_PATH_example.sln
LAI_PATH_example.vcxproj

4. GSL - GNU Scientific Library for windows
gsl\
cblas.lib
cblas32.lib
gsl.lib
gsl32.lib

5. Executable file of Path Length Distribution Model
LAI_PATH.exe		(64bit version)
LAI_PATH_win32.exe	(32bit version)

6. Examples for input files:
input_example1.txt
input_example2.txt
input_example3.txt

Usages:
LAI_PATH -i in.txt -o out.txt
LAI_PATH -i in.txt
LAI_PATH -h

For any question, bug report, please contact Ronghai HU ( huronghai@ucas.edu.cn ).






