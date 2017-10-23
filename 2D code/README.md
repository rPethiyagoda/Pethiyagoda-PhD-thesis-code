
# Pethiyagoda-PhD-thesis-code
Sample code used to generate a figure from:

Pethiyagoda, R., McCue, S. W., & Moroney, T. J. (2017). Efficient computation of two‐dimensional steady free‐surface flows. Int. J. Num. Meth. in Fluids.

------------------------------------------------------------------------
## PROGRAM REQUIREMENTS
------------------------------------------------------------------------

#### Required:


Matlab

Sundials KINSOL (sundialsTB Matlab toolbox in v2.6.2): https://computation.llnl.gov/projects/sundials/sundials-software

#### Optional:

Intel MKL

CUDA 5 compatible GPU

------------------------------------------------------------------------
## Compiling mex files
------------------------------------------------------------------------

If you wish to use the Intel MKL library you must use the mex files. If the included binaries fail, there is an example makefile used to compile the mex files. The directories for MATLAB and the Intel compiler must be changed to reflect your system.

The mex source files are:
- factorBand.c
- factorDense.c
- solveBand.c
- solveDense.c

In order to compile factorBand.c run 'make NAME=factorBand' for example.

------------------------------------------------------------------------
## SCRIPT FILES - run these
------------------------------------------------------------------------

runExample.m - Generates figure 9 of Pethiyagoda et al. (2017)

------------------------------------------------------------------------
## OTHER MAIN FILES
------------------------------------------------------------------------

computeSurface.m     - Computes the solution surfaces

F2D.m	       - The system of nonlinear equations

plotSurf.m	         - Plots a surface

NonlinearBlockJacPre.m  - Create the preconditioner matrix (dense sotrage)

NonlinearBandJacPre.m   - Create the preconditioner matrix (banded sotrage)

------------------------------------------------------------------------
## AUXILIARY FILES
------------------------------------------------------------------------

#### 2014a paper

bMultnonlin.cu	      - Perform the matrix vector product for the B sub matricies of the preconditioner on the GPU

cMultnonlin.cu	      - Perform the matrix vector product for the C sub matricies of the preconditioner on the GPU

factorBand.c	    - Factor a matrix with banded storage (intel MKL, source)

factorBand.mexa64 - Factor a matrix with banded storage (intel MKL, compiled)

factorDense.c	    - Factor a matrix with dense storage (intel MKL, source)

factorDense.mexa64- Factor a matrix with dense storage (intel MKL, compiled)

fGutsGPU.cu	      - Body of the system of nonlinear equations to be evaluated on a GPU

getXY.m	      - Compute the free surface and bottom coordinates

intWeight.m	      - Generate weighting vector for numerical integration

preSolveBlock.m	  - Apply the inverted preconditioner matrix to a vector (Dense sotrage)

preSolveBnd.m	    - Apply the inverted preconditioner matrix to a vector (Banded sotrage)

setTo0.m	        - Set global values to 0

solveBand.c	      - Solve Ax=b with banded storage (intel MKL, source)

solveBand.mexa64  - Solve Ax=b with banded storage (intel MKL, compiled)

solveDense.c	    - Solve Ax=b with dense storage (intel MKL, source)

solveDense.maxa64 - Solve Ax=b with dense storage (intel MKL, compiled)

startup_STB.m	    - Initiate KINSol in MATLAB


