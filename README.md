# ME382 - Final Project
This is the final project of ME382 - Computational Methods in 
Thermo-Fluid Science. It solves 2D axisymmetric cylindrical heat
conduction problem with convective boundary conditions. Uniform heat generation or Arrhenius type heat generation can be added. Anisotropic heat conductivity can be defined. (Uniform kr and uniform kz) Different convective heat transfer coefficients can be defined on all three surfaces. 

This solver is based on implicit finite difference method. Solution of the linear system can be done by both exact solver of numpy, or iterative GMRES solver with iLU preconditioner. Arrhenius heat generation is treated with Picard's iteration. Analytical solution based on seperation of variables is available for no generation case. The grid is structured, the domain divided into equally spaced cells based on number of grids provided for each direction.

Running the code:
The main.py script should be run such as "python main.py" or "python3 main.py". This code is compatible with python3, not python2 (older versions). This code requires the installation of following packages:

* numpy
* scipy
* pandas
* time 
* csv
* matplotlib (for postprocessing only)
* ffmpeg (for contourplot animation movie only)

The code consist of multiple scripts. Explanation of which are provided below:

* main.py
This is the main routine to run. All input parameters, solution and postprocessing options should be manually adjusted here. It evaluates the necessary parameters for the grid and calls the solvers, writes the parameters and calls postprocessing tools if desired. 

* finiteDifference.py
This script is the core finiteDifference solver. It solves the entire domain over time and returns all temperatures and the average temperatures. Two subroutines are available, one is for uniform generation (or no generation) and the other is for Arrhenius type generation.

* analytical_nogen.py
This solves the no generation case with seperation of variables method. It uses 9 terms if Fo<0.2 and 6 terms if Fo>0.2. Current version is accurate, but slow. To speed up, an adaptive selection algorithm can be implemented in the future for decision of number of terms. It returns all temperatures and average temperatures in same structure as finiteDifference.py gives.

* lumped.py
This solves Arrhenius heat generation case with a lumped model.

* lineplot.py
It contains many subroutines for lineplotting, such as: radial distribution, axial distribution, average temperature over time, surface temperature over time...

* contourplot.py
At any time, a contour plot of analytical solution and finite difference solution can be obtained by using this script.

* createanimation.py
Contour plots of any two solutions can be animated over time in same video for comparison. 


