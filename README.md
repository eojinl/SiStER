# SiStER
This repository includes the code necessary to reproduce the viscoelastic experiments in: 

Lai, C.-Y., Rines, J. H., Coffey, N. B., Lee, E., Clerc, F., Hilditch, J. P., Gagliardini, O., & L. Stevens. (2026). From blistering to sliding: Ice-sheet stress perturbations from water injection to the bed

Modified from: 

Simple Stokes solver with Exotic Rhologies

An easy-to-use MATLAB code to simulate long-term lithosphere and mantle deformation.

J.-A. Olive, B.Z. Klein, E. Mittelstaedt, M. Behn, G. Ito, S. Howell
jaolive <at> ldeo.columbia.edu


To run the numerical experiment, run the function "SiStER_MAIN_viscoelastic.m": 

1. Parameters can be modified in the input file "SiStER_Input_File_viscoelastic.m"
2. In order to run a Newtonian/Non-Newtonian experiment, the "nonNewtonian" variable in the input file can be toggled false/true
3. The pressure and slope are modified using the variables "BC.p" and "epsilon" respectively

To reproduce the figures: 

1. The file "PNAS_visualization" contains all relevant code 
2. To load all necessary variables, it is recommended to run "SiStER_MAIN_viscoelastic.m" for one iteration
3. The first section in the file "PNAS_visualization" should be run before running any of the other code. All other blocks do not have a particular sequence
