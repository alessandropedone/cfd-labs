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

In this lab, we tackle the following problems:
1. Solve the Stokes problem with GMRES and ILU preconditioner, which means adding parameters
    ```python
    parameters = {'ksp_type': 'gmres', 'pc_type': 'ilu', 'ksp_rtol': 1.e-8, 'ksp_max_it': 10000}
    ```
2. Convergence test for Stokes in the case we know the analytical solution (so we take exactly the case of the last part of the previous lab), both for $\mathbb{P}^{k+1}/\mathbb{P}^{k}$ and $\mathbb{P}^1_b/\mathbb{P}^1$

Additional notes:
- To impose $\int _\Omega p_h \approx 0$ we need to inform the solver about the null space for Dirichlet boundary conditions
    ```python
    nullspace = MixedVectorSpaceBasis(
        W, 
        [
            W.sub(0), 
            VectorSpaceBasis(constant=True, comm=COMM_WORLD)
        ])
    solve(a == L, w, bcs=bc, nullspace=nullspace)
    ```
- To get the number of iteration from GMRES, you need to use `LinearVariationalSolver` to solve the problem.
- Decreasing $h$ corresponds to more itereations for GMRES.
- The theoretical orders of convergence are
    $$
    \begin{align*}
    \|u - u_h\|_{L^2} &= 
    \begin{cases}
    O(h^{k+1}) \quad &\text{for } \mathbb{P}^{k}/\mathbb{P}^{k-1}\\
    O(h^{2}) \quad \quad &\text{for }\mathbb{P}^{1}_b/\mathbb{P}^{1}
    \end{cases}\\
    \|\nabla u - \nabla u_h\|_{L^2} &= 
    \begin{cases}
    O(h^{k}) \quad \quad &\text{for } \mathbb{P}^{k}/\mathbb{P}^{k-1}\\
    O(h^{1}) \quad \quad &\text{for }\mathbb{P}^{1}_b/\mathbb{P}^{1}
    \end{cases}\\
    \|p - p_h\|_{L^2} &= 
    \begin{cases}
    O(h^{k}) \quad \; &\text{for } \mathbb{P}^{k}/\mathbb{P}^{k-1}\\
    O(h^{3/2}) \quad \; &\text{for }\mathbb{P}^{1}_b/\mathbb{P}^{1}
    \end{cases}
    \end{align*}
    $$
- We check the orders obtained numerically graphically (log-log plot) and quantiatively (logaritmic comparison)
- We may obtain a worse order for the pressure, since the condition $\int _\Omega p_h \approx 0$ is not precise

bubble convergence

## Lab 4