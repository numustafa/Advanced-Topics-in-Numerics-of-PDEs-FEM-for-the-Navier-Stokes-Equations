# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Purpose

This is a course materials repository for two related TU Berlin mathematics courses:

1. **Numerische Mathematik II fur Ingenieurwissenschaften** (Winter term) — prerequisite course on numerical mathematics, taught by Prof. Dr. Tobias Breiten and Prof. Dr. Sandra May. Covers finite difference methods for model problems (Poisson, heat, linear advection equations).

2. **Advanced Topics in Numerics of PDEs: FEM for the Navier-Stokes Equations** (SS2026) — the main course, taught by Dr. Peter Munch. Focuses on finite element methods applied to increasingly complex PDE problems, culminating in the Navier-Stokes equations.

## Repository Structure

- `num_math/` — Lecture slides (PDFs) and reference images for the prerequisite Numerische Mathematik II course (Lectures 01-20, intro slides)
- `week_01(Motivation_and_FEM)/` — Week 1 materials for the FEM/Navier-Stokes course: Chapter 1 (Motivation), Chapter 2 (Finite Element Computations), and Homework 0

## Course Syllabus (FEM for Navier-Stokes)

The course follows this progression (from Chapter 2 table of contents):
1. Motivation
2. Finite element computations
3. Software development with Linux and C/C++
4. Elliptic & parabolic problems
5. Hyperbolic problems
6. Linear solvers & multigrid
7. Stokes problem
8. Navier-Stokes problem
9. Adaptive mesh refinement
10. High-performance computing

## Software Stack

The course uses:
- **deal.II** (v9.6.2) — C++ finite element library, the primary computational tool
- **CMake** — build system for deal.II projects
- **Gmsh** — computational mesh generation
- **ParaView** — postprocessing/visualization of simulation results
- **meshpy** — Python library for unstructured mesh generation (used in HW_00 for 2D triangular meshes)
- **gcc/clang** — C/C++ compilers
- Reports are written in **LaTeX**

## Build Pattern for deal.II Projects

```bash
cd <project_dir>
mkdir -p build && cd build
cmake ..                    # or cmake -DDEAL_II_DIR="<path>" ..
make
./<executable>
```

## Content Format

All materials are PDFs and images — there is no source code in this repository yet. Future weeks will likely add C++ source code for deal.II-based FEM implementations and Python scripts for mesh generation and analysis.
