# SFF – Interacting SSH Model

This repository provides Julia scripts for computing and analyzing the **Spectral Form Factor (SFF)** in an interacting Su-Schrieffer-Heeger (SSH) model with disorder and sublattice symmetry. Both topological and trivial phases are investigated.

---

## Repository Structure

### 1. `sff_test.jl` – SFF Data Generation

Main script for computing the SFF (sequential code).

**Key Parameters:**
- `L = 20`           # Number of lattice sites
- `Np = 4`            # Number of particles (filling)
- `V = 2.0`            # Interaction strength
- `τ = 1.0`            # Base hopping amplitude
- `Δ = -0.1`           # Dimerization parameter  
- `t1_top = (τ + Δ) / 2.0` # Topological hopping
- `t2_top = (τ - Δ) / 2.0` # Topological hopping
- `NR = 10`            # Realizations per disorder config
- `σ = 0.01`           # Disorder strength
- `N_samples = 100`      # Number of disorder samples
- `time_points = 500`
- `time_array = exp10.(range(-1, 4, length=time_points))`
- `β = 0.5`            # Inverse temperature

**Experimental Variables:**
- `V` — Interaction potential (vary it slowly near phase transition point)
- `Np` — Particle filling (flexible, can be changed on GPU)
- `Δ` — Test both `+Δ` and `−Δ` values

---

### 2. `sff_plots.jl` – Plotting & Ramp Slope Extraction

- Loads output data from `sff_test.jl`
- Plots SFF as a function of time
- Extracts ramp slope from the linear regime

---

### 3. GPU Code for SFF

Scripts for GPU-accelerated computation. Simply run and test for performance and accuracy of - `sff_gpu_test_1.jl`

---

### 4. Eigenvalue Comparison Scripts (No Randomization in \( w_{ij} \))

Compare eigenvalues obtained from:
- **GPU-accelerated computation** (e.g., NVIDIA Ada 5000)
- **Sequential CPU computation**

**Scripts:**
- `ssh_cuda_eig.jl`
- `ssh_sq_eig.jl`

**Data Files:**
- `eigenvalues_L20_Np4_gpu.txt`
- `eigenvalues_sq_20_4.txt`
- `difference_plot_1.png` (scatter plot of absolute differences between eigenvalues)

---

## Usage

1. **Run `sff_test.jl`** to generate SFF data.
2. **Visualize and analyze** using `sff_plots.jl`.
3. **Compare GPU and CPU eigenvalues** with the provided eigenvalue scripts and plot.

---

## Contact

For questions or collaboration, please contact [Theoreticalmind](https://github.com/Theoreticalmind).
