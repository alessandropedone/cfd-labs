# CFD labs
This repository contains the solutions of the lab sessions of __COMPUTATIONAL FLUID DYNAMICS COURSE (A.Y. 2025/2026)__, held by Prof. Nicola Parolini and Prof. Lorenzo Valdettaro.

In particular, each folder corresponds to a lab and contains the problem, the solution, the reference code provided by the professor (which may not overlap perfectly with the questions of the lab) and my solution.

In this README, you can find a quick overview of the content of each lab.

## Lab 1. Firedrake.
This lab contains a general introduction to Firedrake with some basic notes about solving numerical problems with this library. Contents:
- Define and plot the mesh
- Define function space and trial/test functions
- Define the source term and the bilinear form
- Add the boundary conditions
- Solve the problem
- Plot the numerical solution
- Algebraic formulation of the problem (with some attention to the algebraic formulation of Dirichlet boundary conditions)

> In the weak formulations specify you need to introduce the lift.

## Lab 2. Stokes: system assembly.
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

There is also a section about the null space of the Stokes monolithic sytem in the case of fully Dirichlet boundary conditions. In particular, if we call $\Sigma$ the matrix associated to the problem with boundary conditions, we observe that $\ker(\Sigma) = \mathrm{span} ([0,1])$, i.e. the space of zero velocity and constant pressure.

> Obvisouly translations of the pressure are allowed we have only Dirichlet boundary conditions on the velocity.

> So, when we have Dirichlet boundary conditions, we want to impose $p \in L^2_0$, i.e. $\int_\Omega p = 0$ at the discrete level. And we see this in next lab.

## Lab 3. Stokes: preconditioning and convergence.

In this lab, we tackle the following problems:
1. Solve the Stokes problem with GMRES and ILU preconditioner, which means adding parameters
    ```python
    parameters = {'ksp_type': 'gmres', 'pc_type': 'ilu', 'ksp_rtol': 1.e-8, 'ksp_max_it': 10000}
    ```
2. Convergence test for Stokes in the case we know the analytical solution (so we take exactly the case of the last part of the previous lab), both for $\mathbb{P}^{k+1}/\mathbb{P}^{k}$ and $\mathbb{P}^1_b/\mathbb{P}^1$

> There is also a starting section about the solving linear systems and preconditioning in Firedrake.

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
- Even if we inform the solver, the presence of a nontrivial null space may spoil the order of convergece we observe, sending us back to the first order: for this reason we introduce a new quantity $$\hat p_h = p_h + \frac{1}{|\Omega|}\int_\Omega (p-p_h)$$, since in this case we know the exact solution $p$
- To get the number of iteration from GMRES, you need to use `LinearVariationalProblem` and `LinearVariationalSolver` to separate the creation of the problem and the solution step
- Decreasing $h$ corresponds to more itereations for GMRES if we use a nonoptimal preconditioner
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

## Lab 4. Stokes: optimal preconditioning.

This lab is about the optimal preconditioners for the Stokes system:

$$
P_d = \begin{bmatrix}
A & 0 \\
0 & \frac{1}{\nu}M_p
\end{bmatrix}, \quad \quad 
P = \begin{bmatrix}
A & 0 \\
\pm B & \frac{1}{\nu}M_p
\end{bmatrix}.
$$

Additional notes:
- We define the preconditioner passing the corresponding bilinear form to `LinearVariationalProblem` as an input arguement
- We use `fieldsplit` to activate block splitting (additive for block diagonal preconditioner and multiplicative for block lower triangular preconditioner), and for each block of the splitting we use different preconditioners (LU or ILU)
- The number of iterations of GMRES is almost constant with respect to $h$ with optimal preconditioners, but when we introduce ILU for one of the blocks we may lose this property

## Lab 5. Stokes: pressure matrix (Schur) method.

The main topic of this lab is the pressure matrix method for Stokes, that corresponds to solving for the Schur complement. Contents:
- Idea: we don't compute the inverse of $A$, but we define a solver with specific PETSc parameters (both fro $A$ and $S$)
- Application: we use this approach to sove the first time-dependent problem for the Stokes equations

Annotations:
- To the define the solver for the Schur complement we rely on 
    ```python
    S = PETSc.Mat().createSchurComplement(A, A, Bt, asB, None)
    ```
- There are several ways to use the pressure mass matrix to precondition the Schur complement (full matrix, the diagonal or the lumped mass matrix)
- We impose weakly Dirichlet boundary conditions, that means using Neumann-type structure
    $$
    u = g + \varepsilon (\nabla u - p)n \implies (\nabla u - p)n  = \frac{1}{\varepsilon} (u-g).
    $$

## Lab 6. NS: Euler's and BDF methods.

This lab contains linear methods to solver time-dependent Navier-Stokes equations. Contents:
- Euler methods: fully explicit, viscous implicit and advection semi-implicit
- BDF methods: fully explicit, viscous implicit and advection semi-implicit
- Convergence tests in time: Euler is first order and BDF is second order

Annotations:
- We can change the preconditioner depeding on the explicit/implicit strategy we choose, since the top left block of the monolithic system changes (see first block of code)
- To observe che convergence for the pressure we need again a correction for its mean

# Lab 7. NS: Picard's and Newton's methods.

# Lab 8. NS: Scalable Nonlinear Equations Solver (SNES).

# Lab 9. NS: Stabilization techniques.

Contents:
- Conditional definition of a parameter:
    ```python
    deltaK = delta * conditional(lt(ReK, one), h/(ubar+1e-10), h**2 * Re)
    ```

# Lab 10. NS: Time dependent simulation with Chorin-Temam method.

The goal of this lab is to implement the incremental Chorin-Temam method to solve the time-dependent Navier-Stokes equation for a flow past a cylinder. Contents:
- Forms for the three steps of the incremental Chorin-Temam method
- Boundary conditions: no-slip on the cylinder, no-slip or free-slip on the walls, Dirichlet condition at the inflow, Neumann condition on the stress at the outflow 
- The initial condition is given by the Stokes equations
- At each time step the drag force and the lift force are computed

# Lab 11. NS: Boussinesq buoyancy model as a coupled problem.

This lab contains the solution of the thermal boundary layer near a hot plate, modeled by the steady incompressible Navier–Stokes equations with the Boussinesq buoyancy model. Contents:
- Fixed point method
- SUPG stabilization for the advection-diffusion of the temperature

# Lab 12. Level-set method for multi-phase flows.

Dir + Free-slip => I still need the null space

# Lab 13. Inexact factorization.

