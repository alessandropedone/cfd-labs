# CFD labs
This repository contains the solutions of the lab sessions of COMPUTATIONAL FLUID DYNAMICS COURSE A.Y. 2025/2026.
In particular, in each folder there are the problem, the solution, the proposed implemented solution by the teachers and my solution (in TEMPLATE files).

In the following there is quick overview of the content of each lab.

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
- Tests with non inf-sup stable spaces  $\mathbf{P}^1/\mathbf{P}^0$ and $\mathbf{P}^1/\mathbf{P}^1$ 
- Tests with inf-sup stable spaaces $\mathbf{P}^2/\mathbf{P}^1$ and $\mathbf{P}_b^1/\mathbf{P}^1$
- Quiver plot for the velocity field
- Compute the $H^1$, $H^1_0$ and $L^2$ errors
- Compute the number of DOF
- Additional parameters for the solver
- Inspection of the matrix (sparsity pattern and difference when boundary conditions are on)

There is also a section about the null space of the Stokes monolythic sytem in the case of fully Dirichlet boundary conditions. In particular, we observe that $\ker(\Sigma) = \mathrm{span} ([0,1])$, i.e. the space of zero velocity and constant pressure.

> Obvisouly translations of the pressure are allowed we have only Dirichlet boundary conditions on the velocity.


## Lab 3
impose L20 for the discrete pressure