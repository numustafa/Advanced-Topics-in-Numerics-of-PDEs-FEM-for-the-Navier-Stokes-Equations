# Advanced Topics in Numerics of PDEs: FEM for the Navier-Stokes Equations

> **SS 2026 | TU Berlin | Institute of Mathematics**
> Taught by Dr. Peter Munch

A comprehensive learning repository for the graduate-level course on Finite Element Methods applied to the Navier-Stokes equations. This repo tracks lecture materials, homework solutions, personal study notes, and (eventually) C++ implementations using the [deal.II](https://www.dealii.org/) finite element library.

---

## Course Overview

The course builds from first principles of FEM to a full Navier-Stokes solver, covering both mathematical theory and hands-on implementation:

| # | Topic | Status |
|---|-------|--------|
| 1 | Motivation & flow problems | Done |
| 2 | Finite element computations | In progress |
| 3 | Software development with Linux and C/C++ | Upcoming |
| 4 | Elliptic & parabolic problems | Upcoming |
| 5 | Hyperbolic problems | Upcoming |
| 6 | Linear solvers & multigrid | Upcoming |
| 7 | Stokes problem | Upcoming |
| 8 | Navier-Stokes problem | Upcoming |
| 9 | Adaptive mesh refinement | Upcoming |
| 10 | High-performance computing | Upcoming |

## Repository Structure

```
.
├── num_math/                          # Prerequisite: Numerische Mathematik II
│   ├── intro_slides.pdf               #   Course introduction (Breiten & May)
│   ├── Lecture01_Oct14.pdf            #   Finite difference methods
│   ├── ...                            #   Lectures 02-20 (FD, heat eq, FEM)
│   ├── Lecture15.pdf                  #   Variational formulation & finite elements
│   ├── Lecture16.pdf                  #   1D/2D element catalog, stiffness matrix
│   ├── Lecture17.pdf                  #   Lax-Milgram, FE analysis
│   └── FEM_topics.png                 #   Key FEM concepts overview
│
├── week_01(Motivation_and_FEM)/       # Week 1: Chapters 1-2
│   ├── chapter_01.pdf                 #   Ch.1: Motivation (flow problems, NS eqs)
│   ├── chapter_02_empty.pdf           #   Ch.2: Finite element computations (fill-in)
│   ├── HW_00.pdf                      #   Homework 0: FEM review + software setup
│   ├── FEM_study_guide.tex            #   Personal study guide (LaTeX source)
│   └── FEM_study_guide.pdf            #   Compiled study guide (16 pages)
│
├── CLAUDE.md                          # AI assistant context file
└── README.md
```

> **Planned directories** (will be added as the semester progresses):
> `week_02/`, `week_03/`, ..., and a `code/` directory for deal.II implementations.

## Prerequisites

This course assumes familiarity with:

- **Weak formulations** and the variational framework for PDEs
- **Basic FEM**: shape functions, assembly, stiffness matrix
- **Galerkin methods**: orthogonality, best approximation, Lax-Milgram
- **Linear algebra**: sparse matrices, iterative solvers

The `num_math/` directory contains the full prerequisite lecture series (Numerische Mathematik II, Breiten & May, WS 2025/26) covering finite difference methods through finite elements from scratch -- Lectures 12-17 are particularly relevant.

## Software Stack

| Tool | Purpose |
|------|---------|
| [deal.II](https://www.dealii.org/) v9.6.2 | C++ finite element library |
| CMake | Build system for deal.II projects |
| gcc / clang | C/C++ compilers |
| [Gmsh](https://gmsh.info/) | Computational mesh generation |
| [ParaView](https://www.paraview.org/) | Postprocessing & visualization |
| LaTeX | Homework reports and study notes |
| Python + meshpy | Mesh generation for introductory exercises |

### Building a deal.II project

```bash
cd <project_dir>
mkdir -p build && cd build
cmake ..                                    # or: cmake -DDEAL_II_DIR="<path>" ..
make
./<executable>
```

### Installing deal.II (Ubuntu / WSL)

**Via package manager (quick):**
```bash
export REPO=ppa:ginggs/deal.ii-9.6.0-backports
sudo apt-get update && apt-get install -y software-properties-common
sudo add-apt-repository $REPO
sudo apt-get update
sudo apt-get install libdeal.ii-dev
```

**From source via candi (full control):**
```bash
git clone https://github.com/dealii/candi.git && cd candi
sed -i 's/^#\?DEAL_II_VERSION=.*/DEAL_II_VERSION=v9.6.2/' candi.cfg
./candi.sh --prefix="$HOME/software/dealii-9.6/" \
           --packages="p4est trilinos dealii" \
           --jobs=4
```

## Study Guide

The `FEM_study_guide.pdf` in `week_01/` is a 16-page personal reference that:

- Walks through the FEM pipeline: **strong form -> weak form -> Galerkin -> linear system**
- Uses notation consistent with the prerequisite lectures ($\psi_i$, $V_n$, $A_n$, $d_i$) with a translation table for the Chapter 2 slides ($\varphi_j$, $V_h$, $K$, $U_j$)
- Covers variational formulations, Galerkin orthogonality, Lax-Milgram, element types, assembly, error estimates, and adaptive mesh refinement
- Includes a prioritized reading roadmap through the prerequisite lectures

To recompile after edits:
```bash
cd "week_01(Motivation_and_FEM)"
pdflatex FEM_study_guide.tex && pdflatex FEM_study_guide.tex
```

## Homework

| HW | Topic | Deadline | Key tasks |
|----|-------|----------|-----------|
| 0 | Setup & review | Sun, May 3 2026 | FEM code review, install terminal/gcc/deal.II |
| 1 | Variational forms | TBD | Derive variational forms, compute mass matrix, extend Python code to quads |

## Acknowledgments

- **Course instructor:** Dr. Peter Munch, Institute of Mathematics, TU Berlin
- **Prerequisite lectures:** Prof. Dr. Tobias Breiten & Prof. Dr. Sandra May, TU Berlin
- **deal.II library:** W. Bangerth, T. Heister, et al. -- [dealii.org](https://www.dealii.org/)
