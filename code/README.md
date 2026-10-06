# Codebase Documentation

This document categorizes and explains the purpose of all scripts and modules in the `code/` directory to help distinguish between core thesis code, active experiments, and deprecated/scratchpad scripts.

---

## 1. Top-Level Scripts (`code/`)
- **`run_all_experiments.m`**: Master script to run FDTD vs PSTD error comparisons.
- **`run_exp1_fdtd.m`, `run_exp1_pstd.m`**: Individual experiment runners for specific wave configurations.
- **`generate_table.m`, `print_table.m`, `table_data.mat`**: Utility scripts to compute and print the error scaling tables across different grid resolutions $N$.
- **`test_size.m`**: Scratchpad for testing array/matrix sizes. *(Candidate for deletion)*

---

## 2. PSTD (Pseudo-Spectral Time-Domain)
Located in `code/pstd/`.

### Core Source (`code/pstd/src/`)
- **`solver.m`**: The baseline monolithic 1D PSTD solver.
- **`ard_solver.m`**: The primary Adaptive Rectangular Decomposition (ARD) solver using the Exact Exponential Strang Splitting (used for the thesis results).
- **`solver_wrapper_ard.m`**: Wrapper aligning ARD with the monolithic interface for easy comparisons.
- **`get_residue_matrix.m`, `modal_propagators.m`**: Matrix assemblers for the PSTD spatial derivatives and exact temporal evolution.
- **`compute_cfl.m`, `compute_parabolic_limit.m`**: Computes the stability limit $\Delta t_{\text{max}}$ for the viscoelastic wave equation.
- **`cosine_interpolation.m`**: Helper for generating smooth initial conditions.

**Deprecated integration operators *(Candidates for deletion)***:
We explored several splitting operators before settling on the Exact Exponential matrix. These have large truncation errors and were abandoned:
- `ard_solver_cn_kick.m`, `ard_solver_without_nu.m`, `ard_solver_no_splitting.m`
- `compute_cfl_cn_kick.m`, `compute_cfl_exp_euler.m`

### Scripts & Experiments (`code/pstd/scripts/`)
- **`regenerate_all_thesis_images.m`**: Master script that runs all 6 experimental wave scenarios, comparing Monolithic vs ARD PSTD, and exports snapshots for the LaTeX document.
- **`run_solver.m`, `run_ard_solver.m`**: Simple runner scripts to visualize a single simulation.
- **`convergence_test.m`, `convergence_test_ard.m`, `convergence_test_ard_space_order.m`**: Measure convergence rates ($\Delta t \to 0$ and $\Delta x \to 0$) for the ARD solver against monolithic reference.
- **`compare_strang_vs_exponential.m`, `compare_strang_euler_vs_strang_cn.m`, `compare_nu_correction.m`**: Scripts used to measure and prove that the Crank-Nicolson/Euler operators introduce massive $O(\Delta t)$ splitting errors, justifying the exponential approach.

**Debugging Scratchpads *(Candidates for deletion)***:
- `debug_cfl.m`, `verify_cfl_empirical.m`, `check_nx_trend.m`, `test_G_matrix.m`, `convergence_test_cn_kick.m`, `convergence_test_exp_euler.m`

---

## 3. FDTD (Finite-Difference Time-Domain)
Located in `code/fdtd/`.

### Core Source (`code/fdtd/src/`)
- **`solver.m`**: The baseline monolithic 1D FDTD solver.
- **`swr_solver.m`**: Schwarz Waveform Relaxation (OSWR) 1D FDTD solver.
- **`get_optimal_robin_params.m`**: Solves the theoretical min-max optimization problem to find the optimal $(q, r)$ Robin transmission parameters.
- **`interior_assembly.m`, `apply_bc.m`**: Matrix assemblers for FDTD differences and boundary conditions.

### Scripts & Experiments (`code/fdtd/scripts/`)
- **`run_solver.m`, `run_swr_solver.m`**: Basic visual runners.
- **`convergence_test.m`, `convergence_test_swr.m`**: Verifies error scaling for the OSWR FDTD solver.
- **`convergence_test_swr_part2.m`**: Plots the min-max convergence factor $\rho$ over frequencies $\omega$.
- **`convergence_test_swr_part3.m`**: Runs the 2D grid sweep of $(q, r)$ to plot the empirical error basin vs theoretical optimum.
- **`convergence_test_swr_part4.m`**: Sweeps physical parameters $(\gamma, \nu)$ to compute a 2D map of the optimal $\rho_{\text{opt}}$.
- **`convergence_test_swr_part5.m`**: Plots parametric continuous trajectories of optimal $(q, r)$ as $(\gamma, \nu)$ vary.
- **`generate_sweeps.m`**: Helper to generate sweeps.

**Debugging Scratchpads *(Candidates for deletion)***:
- `code/fdtd/test_1d.m`

---

## 4. FDTD 2D (Sandbox)
Located in `code/fdtd_2d/`.

Contains the 2D extension of the FDTD OSWR solver.
- **`src/`**: `solver_2d.m`, `swr_solver_2d.m`, `interior_assembly_2d.m`, `apply_bc_2d.m`.
- **`convergence_test_swr_2d.m`, `generate_error_plot_2d.m`**: Scripts running the 2D min-max experiments.

**Debugging Scratchpads *(Candidates for deletion)***:
- `test_err.m`, `test_fixed.m`, `test_initial.m`, `test_interp.m`, `fdtd_2d_test.m`, `fdtd_2d_test_damped.m`.

---

## 5. Differentiable Solvers (DPNN)
Located in `code/dpnn/`.
- **`dpnn_1d_compare_solvers.ipynb`**: PyTorch Jupyter notebook training the SIREN network to reconstruct initial conditions using Differentiable FDTD and PSTD solvers.
- **`diff_check.py`**: A python script verifying the inverse parameter estimation gradients.
- **`run_notebook.bat`**: Simple launch script.
