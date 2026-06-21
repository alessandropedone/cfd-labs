# CFD labs
This repository contains the solutions of the lab sessions of __COMPUTATIONAL FLUID DYNAMICS COURSE A.Y. 2025/2026__.

In particular, each folder corresponds to a lab and contains the problem, the solution, the reference code provided by the professor (which may not overlap perfectly with the questions of the lab) and my solution.

In this README, you can find a quick overview of the content of each lab.

## Lab 1
This lab contains a general introduction to Firedrake with some basic notes about solving numerical problems with this library. Contents:
- Define and plot the mesh
- Define function space and trial/test functions
- Define the source term and the bilinear form
- Add the boundary conditions
- Solve the problem
- Plot the numerical solution
- Algebraic formulation of the problem (with some attention to the algebraic formulation of Dirichlet boundary conditions)

> In the weak formulations specify you need to introduce the lift.

## Lab 2
This lab contains the Stokes system assembly and solution, in very simple cases (Coutte flow and Poiseuille flow). Contents:
- Define and plot a rectangular mesh
- Define vector and mixed function space
- Tests with non inf-sup stable spaces  $\mathbb{P}^1/\mathbb{P}^0$ and $\mathbb{P}^1/\mathbb{P}^1$ 
- Tests with inf-sup stable spaaces $\mathbb{P}^2/\mathbb{P}^1$ and $\mathbb{P}_b^1/\mathbb{P}^1$
- Quiver plot for the velocity field
- Compute the $H^1$, $H^1_0$ and $L^2$ errors
- Compute the number of DOF
- Additional parameters for the solver
- Inspection of the matrix (sparsity pattern and difference when boundary conditions are on)

There is also a section about the null space of the Stokes monolythic sytem in the case of fully Dirichlet boundary conditions. In particular, if we call $\Sigma$ the matrix associated to the problem with boundary conditions, we observe that $\ker(\Sigma) = \mathrm{span} ([0,1])$, i.e. the space of zero velocity and constant pressure.

> Obvisouly translations of the pressure are allowed we have only Dirichlet boundary conditions on the velocity.

> So, when we have Dirichlet boundary conditions, we want to impose $p \in L^2_0$, i.e. $\int_\Omega p = 0$ at the discrete level. And we see this in next lab.

## Lab 3

impose L20 for the discrete pressure: inform the solver with the null space 

First exercise to do (cavity problem)

LinearVariationalSolver to get number of iterations

Parameters for GMRES with tolerance $1e-8$ with ILU preconditioner
```python
parameters = {'ksp_type': 'gmres', 'pc_type': 'ilu', 'ksp_rtol': 1.e-8, 'ksp_max_it': 10000}
```

Decrease h implies more itereations for GMRES

bubble convergence

convergence test with fully Dirichlet conditions (same problem as Lab 2) and theorical orders of convergence for pressure and velocity review

plots and computation of oder for convergence test

## Lab 4