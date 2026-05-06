# CFD-Julia: 12 Steps to Navier-Stokes in Julia

[![Julia](https://img.shields.io/badge/Julia-1.x-9558B2?style=flat&logo=julia&logoColor=white)](https://julialang.org/)
[![Jupyter](https://img.shields.io/badge/Jupyter-Notebooks-F37626?style=flat&logo=jupyter&logoColor=white)](https://jupyter.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A Julia implementation of the classic **"12 Steps to Navier-Stokes"** curriculum, originally developed by [Prof. Lorena Barba's group](https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/) in Python. This project was developed as part of an academic course/thesis in Computational Fluid Dynamics (CFD), progressively building the numerical tools needed to solve the incompressible Navier-Stokes equations from scratch.

---

## Table of Contents

- [Overview](#overview)
- [Steps](#steps)
- [Extended Features](#extended-features)
- [Requirements](#requirements)
- [Getting Started](#getting-started)
- [Project Structure](#project-structure)
- [Acknowledgements](#acknowledgements)

---

## Overview

This repository walks through 12 incremental steps — from simple 1D linear convection all the way to full 2D Navier-Stokes cavity and channel flow — implemented entirely in **Julia**. Each step builds on the last, introducing new physics and numerical methods along the way.

Beyond the core 12 steps, the project includes:

- **3D Navier-Stokes** simulation extending the 2D formulation to three dimensions
- **GPU acceleration** (CUDA via `CUDA.jl`) for Step 1, with a benchmarking comparison
- **3D Channel Flow** as a standalone simulation
- **Visualizations** built with [Makie.jl](https://makie.juliaplots.org/)

---

## Steps

| Step | Topic | Dimension |
|------|-------|-----------|
| 1  | Linear Convection | 1D |
| 2  | Nonlinear Convection | 1D |
| 3  | Diffusion | 1D |
| 4  | Burgers' Equation | 1D |
| 5  | Linear Convection | 2D |
| 6  | Nonlinear Convection | 2D |
| 7  | Diffusion | 2D |
| 8  | Burgers' Equation | 2D |
| 9  | Laplace Equation | 2D |
| 10 | Poisson Equation | 2D |
| 11 | Cavity Flow (Navier-Stokes) | 2D |
| 12 | Channel Flow (Navier-Stokes) | 2D |

Each step is organized in its own folder and implemented as a Julia script (`.jl`) and/or Jupyter Notebook (`.ipynb`) with inline comments explaining the numerical methods.

---

## Extended Features

### 3D Navier-Stokes (`3DNavier.jl`)
Extends the 2D formulation to a full three-dimensional simulation of incompressible flow, solving the Navier-Stokes equations on a 3D domain.

### 3D Channel Flow (`Channel Flow 3D/`)
A dedicated simulation of pressure-driven channel flow in 3D, serving as a more realistic benchmark case beyond the classical 2D formulation.

### GPU Acceleration (`GPU Step 1/`, `GPUKernel.jl`)
Step 1 (Linear Convection) has been ported to run on the GPU using `CUDA.jl` with custom CUDA kernels. Use `Benchmark.jl` to compare CPU vs. GPU performance:

```julia
julia Benchmark.jl
```

### Visualizations with Makie (`PlotMakie.jl`, `PlotNavier.jl`)
Results are visualized using `Makie.jl` / `CairoMakie.jl`, producing publication-quality plots of velocity fields, pressure distributions, and flow streamlines.

---

## Requirements

- [Julia 1.x](https://julialang.org/downloads/)
- [IJulia.jl](https://github.com/JuliaLang/IJulia.jl) — to run Jupyter notebooks
- [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) — required only for GPU steps (needs an NVIDIA GPU)
- [Makie.jl](https://github.com/MakieOrg/Makie.jl) / `CairoMakie.jl` — for visualizations
- [Plots.jl](https://github.com/JuliaPlots/Plots.jl) — used in some steps

Install Julia dependencies from the Julia REPL:

```julia
using Pkg
Pkg.add(["IJulia", "CUDA", "CairoMakie", "Plots"])
```

---

## Getting Started

1. **Clone the repository**

```bash
git clone https://github.com/JesusNagao/CFD-Julia.git
cd CFD-Julia
```

2. **Run a step as a Julia script**

```bash
julia "Step 1 Linear Convection/step1.jl"
```

3. **Open a Jupyter Notebook**

```bash
jupyter notebook
```

Then navigate to the desired step folder and open the `.ipynb` file.

4. **Run the GPU benchmark** (requires NVIDIA GPU)

```bash
julia Benchmark.jl
```

---

## Project Structure

```
CFD-Julia/
├── Step 1 Linear Convection/
├── Step 2 Nonlinear Convection/
├── Step 3 Diffusion/
├── Step 4 Burgers' Equation/
├── Step 5 2D Linear Convection/
├── Step 6 2D Convection/
├── Step 7 2D Diffusion/
├── Step 8 2D Burger's Equation/
├── Step 9 2D Laplace Equation/
├── Step 10 2D Poisson Equation/
├── Step 11 Cavity Flow/
├── Step 12 Channel Flow/
├── Channel Flow 3D/
├── GPU Step 1/
├── Box/
├── 3DNavier.jl          # 3D Navier-Stokes solver
├── Benchmark.jl         # CPU vs GPU benchmark
├── GPUKernel.jl         # Custom CUDA kernels
├── PlotMakie.jl         # Makie visualizations
└── PlotNavier.jl        # Navier-Stokes flow plots
```

---

## Acknowledgements

This project is inspired by and based on the **CFD Python: 12 Steps to Navier-Stokes** course by [Prof. Lorena A. Barba](https://lorenabarba.com/blog/cfd-python-12-steps-to-navier-stokes/). The original course uses Python and NumPy; this repository reimplements the curriculum in Julia, taking advantage of Julia's performance and scientific computing ecosystem.
