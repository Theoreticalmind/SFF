# SFF - Interacting SSH Model

This repository contains Julia scripts to compute and analyze the Spectral Form Factor (SFF) for an interacting SSH model with disorder and sublattice symmetry.We study both topological and trivial phases by varying the interaction strength 𝑉, filling 𝑁𝑝 and dimerization parameter Δ(positive & negative values with the same magnitude).


1.sff_test.jl - sff data generation
 This is the main script for computing the sff . It's a sequantial code.
 Key Paramters:
L = 20          # Number of lattice sites
Np = 4          # Number of particles (filling)
V = 2.0         # Interaction strength
τ = 1.0         # Base hopping amplitude
Δ = -0.1        # Dimerization parameter

t1_top = (τ + Δ) / 2.0   # Topological hopping t1
t2_top = (τ - Δ) / 2.0   # Topological hopping t2

NR = 10                  # Number of realizations per disorder config
σ = 0.01                 # Disorder strength
N_samples = 100          # Number of disorder samples
time_points = 500
time_array = exp10.(range(-1, 4, length=time_points))
β = 0.5                  # Inverse temperature

What to vary for experiments:
. V — interaction potential
. Np — particle filling (can be changed flexibly on GPU)
. Δ test both +Δ and −Δ

2.sff_plots.jl - Plotting and Ramp Slope Extraction
. Reads the output data from sff_test.jl
. Plots SFF vs time
. Extracts ramp slope from the linear regime


This is the sequential code i generated and 

3. GPU code for SFF - you need to just test












4. Eigenvalue Comparison Scripts( without randomization w_ij)
These scripts compare eigenvalues obtained from:
. GPU-accelerated computation (e.g., NVIDIA Ada 5000)
. Sequential CPU computation

These scripts are:







