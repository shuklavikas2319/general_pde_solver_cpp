# Newton–BTDMA Solver for Coupled 1-D Reaction–Diffusion PDEs

A pure C++ implementation of an **implicit finite-difference solver** for a system of two coupled reaction–diffusion equations.

The nonlinear system that arises at each time step is solved by the **Newton–Raphson method**, and the resulting block-tridiagonal linear systems are solved efficiently with a **Block Thomas Algorithm (BTDMA)**.

---

## Governing Equations

The code currently solves the coupled system

$$
\begin{aligned}
\frac{\partial u}{\partial t} &= \frac{\partial^2 u}{\partial x^2} + v, \\
\frac{\partial v}{\partial t} &= \frac{\partial^2 v}{\partial x^2} + u,
\end{aligned}
$$

on the unit interval \(x \in [0,1]\).

**Boundary conditions** (Dirichlet):

* Left (\(x=0\)): \(u=0\), \(v=0\)
* Right (\(x=1\)): \(u=1\), \(v=1\)

**Initial condition**: Zero fields (\(u=v=0\)).

The spatial operators are discretised with second-order central differences on a uniform grid.

Time integration is fully implicit (backward Euler). At each time step, the nonlinear residual is driven to zero by Newton iteration; the Jacobian is obtained by numerical differentiation, and the resulting block-tridiagonal system is solved by BTDMA.

---

## Features

* Fully implicit time marching
* Newton–Raphson nonlinear solver with configurable tolerance and maximum iterations
* Numerically computed Jacobian (easy to replace with analytic derivatives)
* Block-tridiagonal matrix algorithm (BTDMA) for the linear systems
* Clean separation of residual evaluation, Jacobian assembly, linear algebra, and the Newton driver
* Pure C++17, no external libraries required

---

## Repository Layout

```text
├── main.cpp                 # Driver: time loop + Newton loop
├── residual_maker.h         # Residual evaluation declarations
├── residual_maker.cpp       # Residual evaluation (interior + BCs)
├── jacobian_maker.h         # Jacobian assembly declarations
├── jacobian_maker.cpp       # Numerical Jacobian assembly
├── newton_solver.h          # Newton solver declarations
├── newton_solver.cpp        # Newton–Raphson wrapper
├── BTDMA.h                  # Block Thomas algorithm declarations
├── BTDMA.cpp                # Block Thomas algorithm
├── linealgbr.h              # Dense linear-algebra declarations
├── linealgbr.cpp            # Dense linear-algebra helpers
├── build.sh                 # Build script
└── README.md                # Project documentation
```

---

## Building

A simple `g++` command is sufficient:

```bash
g++ -std=c++17 -O2 -o newton_btdma \
    main.cpp residual_maker.cpp jacobian_maker.cpp \
    newton_solver.cpp BTDMA.cpp linealgbr.cpp
```

Alternatively, the source files can be incorporated into a CMake or Makefile project.

---

## Running

On Linux or a compatible shell:

```bash
./newton_btdma
```

On Windows, run the executable produced by your compiler:

```powershell
.\newton_btdma.exe
```

The program prompts for four parameters:

| Prompt                       | Meaning                                    | Typical value |
| ---------------------------- | ------------------------------------------ | ------------- |
| Number of variables          | \(n\) (must be 2 for the current residual) | `2`           |
| Number of nodes              | \(N\) (spatial grid points)                | `51`–`201`    |
| Number of residual equations | Must equal \(n\)                           | `2`           |
| Target time                  | Final simulation time                      | `0.5`–`5.0`   |

Example session:

```text
give the value of number of variables
2
give the value of number of nodes
101
give the value of number of residual equations
2
give the value of target time time_target
1.0
```

During the run, the Newton residual history is displayed for every time step. At the end, the solution vectors \(u\) and \(v\) are printed.

---

## Algorithm Outline

### 1. Residual Evaluation

At every interior node, the discrete residuals of the two PDEs are formed.

Boundary nodes enforce the boundary conditions through a general Robin form:

$$
\alpha\phi+\beta\frac{\partial\phi}{\partial x}=\theta.
$$

### 2. Jacobian Assembly

Each residual equation is differentiated numerically with respect to the unknowns at the left, centre, and right nodes using finite-difference perturbations.

The perturbation size is \(10^{-9}\).

### 3. Newton Step

The linear system

$$
J\,\delta\mathbf{U}=-\mathbf{R}
$$

is solved by BTDMA, and the solution is updated:

$$
\mathbf{U}^{k+1}=\mathbf{U}^{k}+\delta\mathbf{U}.
$$

### 4. Time Advance

After Newton convergence, the old-time level is overwritten, and the next time step begins.

---

## Extending the Code

* **Different physics:** Edit the three residual functions in `residual_maker.cpp` (`residual`, `residual_left`, and `residual_right`).
* **Analytic Jacobian:** Replace the finite-difference loops in `jacobian_maker.cpp`.
* **Boundary conditions:** Change the `BCSpec` arrays used in `main.cpp` and `res_func`.
* **More variables:** The infrastructure supports a configurable number of unknowns \(n\), but the residual definitions must also be generalised beyond the current two-equation system.

---

## Notes and Limitations

* The current residual is hard-coded for two equations; the driver therefore expects `n == num == 2`.
* Matrix inversion inside BTDMA uses classical Gaussian elimination without pivoting. This may be adequate for well-conditioned Jacobians but is not suitable for arbitrary indefinite systems.
* No adaptive time-stepping or spatial adaptivity is implemented.
* Output is currently text-only; CSV or VTK output could be added.

---

## License

No license has been specified yet. Add an appropriate license before distributing the code under explicit open-source terms.

---
