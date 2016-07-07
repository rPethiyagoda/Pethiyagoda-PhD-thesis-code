# Pethiyagoda-PhD-thesis-code
Sample code used to generate figures from:

Pethiyagoda, R., McCue, S. W., Moroney, T. J., & Back, J. M. (2014a). Jacobian-free Newton-Krylov methods with GPU acceleration for computing nonlinear ship wave patterns. J. Comput. Phys., 269, 297--313.

Pethiyagoda, R., McCue, S. W., & Moroney, T. J. (2014b). What is the apparent angle of a Kelvin ship wave pattern? J. Fluid Mech., 758, 468--485.

and

Pethiyagoda, R., McCue, S. W., & Moroney, T. J. (2016). Spectrograms of ship wakes: identifying linear and nonlinear wave signals. Submitted to J. Fluid Mech.

------------------------------------------------------------------------
PROGRAM REQUIREMENTS
------------------------------------------------------------------------

Required:
--------------------------------

Matlab

Sundials KINSOL (sundialsTB Matlab toolbox): https://computation.llnl.gov/casc/sundials/main.html

Optional:
------------------------------------------------------------------------

Intel MKL

CUDA 5 compatible GPU

------------------------------------------------------------------------
Compiling mex files
------------------------------------------------------------------------

If the included binaries fail, there is an example makefile used to compile the mex files.

The directories for MATLAB and the Intel compiler must be chagned to reflect your system.

The mex source files are:

factorBand.c

factorDense.c

solveBand.c

solveDense.c

In order to compile factorBand.c run 'make NAME=factorBand' for example.

------------------------------------------------------------------------
SCRIPT FILES - run these
------------------------------------------------------------------------

runJCP.m 	  - Generates figure 5(b) of Pethiyagoda et al. (2014a)

runJFMwakeAngle.m - calculates the apparent wake angle for a flow pasta dipole with F=0.9 (part of figure 7(b) from Pethiyagoda et al. (2014b))

runJFMDipole.m 	  - Generates solutions to be used by runJFMwakeAngle.m

------------------------------------------------------------------------
OTHER MAIN FILES
------------------------------------------------------------------------

computeSurface.m     - Computes the solution surfaces

BIFunction.m	       - The system of nonlinear equations

plotSurf.m	         - Plots a surface

measureAllPeaks.m    - Measure the apparent wake angle of a series solutions

getPeakAngle.m	     - Measure the apparent wake angle of a single solution

LinearBlockJacPre.m  - Create the preconditioner matrix (dense sotrage)

LinearBandJacPre.m   - Create the preconditioner matrix (banded sotrage)

------------------------------------------------------------------------
AUXILIARY FILES
------------------------------------------------------------------------

baseMult.cu	      - Perform the matrix vector product for the B/C sub matricies of the preconditioner on the GPU

baseMultCPU.m	    - Perform the matrix vector product for the B/C sub matricies of the preconditioner on the CPU

calcMem.m	        - Calculate the memory required for the preconditioner matrix

factorBand.c	    - Factor a matrix with banded storage (intel MKL, source)

factorBand.mexa64 - Factor a matrix with banded storage (intel MKL, compiled)

factorDense.c	    - Factor a matrix with dense storage (intel MKL, source)

factorDense.mexa64- Factor a matrix with dense storage (intel MKL, compiled)

fGutsGPU.cu	      - Body of the system of nonlinear equations to be evaluated on a GPU

getValues.m	      - Compute zeta, phi, and their dervivatives from the vector of unknowns

intWeight.m	      - Generate weighting vector for numerical integration

preSolveBlock.m	  - Apply the inverted preconditioner matrix to a vector (Dense sotrage)

preSolveBnd.m	    - Apply the inverted preconditioner matrix to a vector (Banded sotrage)

setTo0.m	        - Set global values to 0

solveBand.c	      - Solve Ax=b with banded storage (intel MKL, source)

solveBand.mexa64  - Solve Ax=b with banded storage (intel MKL, compiled)

solveDense.c	    - Solve Ax=b with dense storage (intel MKL, source)

solveDense.maxa64 - Solve Ax=b with dense storage (intel MKL, compiled)

startup_STB.m	    - Initiate KINSol in MATLAB



