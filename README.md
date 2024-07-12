Abstract – The study of Cahn-Hilliard equation, a fundamental model that describes the phase separation process in multi-component mixtures. We are using the series and parallel methods to find the best solution. Demonstrating the solution of an equation till its fifth order. The Non-linear / serial method was solved using the required formula whereas for the Parallel method we used certain processors, which uses the split and combine method. We noted that using the parallel method was a better approach and showed better results. Analysis was performed using the open-source package FEniCSx on Python. The program was run using a preconditioned Krylov subspace method for higher-order function spaces. The Krylov subspace solver drastically reduces computational time. The time taken for the execution of each order was recorded and presented. 

 

Index Terms – Cahn-Hilliard; Series Solver; Parallel Solver; FEniCSx; Partial Differential equations 

 

I.  Introduction 

The Cahn-Hilliard equation is a parabolic equation and is typically used to model phase separation in binary mixtures. It involves first-order time derivatives, and second- and fourth-order spatial derivatives. Its physical applications have been extended to many areas of scientific fields such as spinodal decomposition, diblock copolymer, image inpainting, multiphase fluid flows, microstructures with elastic inhomogeneity, tumor growth simulation, and topology optimization [1, 2]. Cahn-Hilliard equation is useful in many applications because it can model the spontaneous separation of a homogeneous mixture into distinct phases over time, it can also capture the evolution of interfaces between different phases and the resulting microstructures.  

Finite element analysis reduces the number of physical tests, lowers the cost, and reduces the amount of time needed for prototyping and testing, all of which speed up design and development [3]. The finite element method (FEM) is a numerical method that finds solutions to field problems that occur in the physical world and several engineering applications. FEM can analyze various issues like fluid dynamics, heat and mass transfer, etc. [4]. The solution obtained is the distribution of the dependent variables over the spatial domain in which the problem is analyzed. These field problems are described using partial differential equations (PDEs) [5]. This work does not focus on the method of the analysis but rather on reducing the computational time required to perform such an analysis. The reduction in computation time is achieved by using parallel processing and Krylov subspace methods to perform the analysis. 

FEniCSx is an open-source program that transforms scientific models into finite element code and applies the FEM to solve PDEs. Although FEniCSx is user-friendly for beginners, it also has significant features for more experienced programmers. The original FEniCS library has been greatly improved and is now known as FEniCSx. Memory parallelization and library improvements are just two of FEniCSx’s many advancements over FEniCS [6, 7]. UFL, Basix, FFCx, and DOLFINx are the four libraries that make up FEniCSx. The computational environment of FEniCSx is DOLFINx.  FEniCSx requires the PDE to be implemented in the Python code in the variational or weak form [8]. Using the mathematical operators of the UFL library, the variational form is written into the Python program. Due to the characteristics mentioned earlier, FEniCSx is quicker, more practical, and able to reach the full range of FEM while yet keeping a solid open-source structure. 

The miniaturization of semiconductor transistors over the last 5 decades has facilitated a massive increase in computational power. This enables the use of multiple processors working in parallel to provide increased processing power. Problems can now be divided into smaller tasks that are computed simultaneously. This reduces computational time and can be used to solve more intensive problems that have not been solved yet. Most commercial FEM softwares rely on robust direct solvers, which makes it difficult to parallelize. Dividing the task across too many processors will be counterproductive and will increase the computational time [9]. For this reason, this work has been restricted to using four processors. 

    Standard Gaussian elimination-based numerical techniques for linear problems are the foundation of direct solvers. These are quite effective in resolving many simpler problems and are advised for problems with up to a few thousand unknowns. Even on the fastest supercomputers, these strategies are incredibly wasteful for very large systems. ILU-preconditioned Krylov subspace solvers solved this issue [10, 11]. Scaling finite element analysis heavily relies on the use of iterative solvers like Krylov solvers. The Krylov solver aims to combine all of the approximations computed up to this point in the iteration process into a better solution [9]. As a result, computing optimum solutions in the Krylov subspace was greatly facilitated. Even when the matrix is very complicated to solve, the Krylov subspace solver achieves the same answer significantly more quickly than the direct solver.  

The iterative Krylov subspace solver also has the added benefit of using far less memory than a direct solver. Compared to solvers used in most FEM software available in the market, this methodology is far more effective due to the improved speed and lower memory usage of Krylov solvers running in parallel across multiple processors. One way of producing more accurate results is to use higher-order function spaces. In this work, we solved PDEs till order 3. However, doing so drastically increases the number of DOFs and hence the computational power and time required to solve the problem. Parallel processing and the Krylov solvers can reduce the time needed to solve such problems. The reduction in computation time is achieved by using parallel processing and Cahn-Hilliard equation to perform the analysis.  The time taken to execute the program under different solver configurations are presented in this paper. 

 

 

II.  EQUATIONS 

A.	Problem definition 

 

The Cahn-Hilliard equation, a parabolic equation typically used in model phase separation in binary mixtures involving first-order time derivatives, and second- and fourth-order spatial derivatives. 

Equation being: 

 

 

 

unknown field is defined as c, the non-convex in c (a fourth-order polynomial is commonly used) is defined as the function f, outward directed boundary normal is defined as n, and M is a scalar parameter. 

B.	Mixed form  

 

Fourth-order equation of Cahn - Hilliard, casting it in a weak form would result in the presence of second-order spatial derivative, and the problem could not be solved using a standard Lagrange finite element basis. One of the solutions is to rephrase the problem as two coupled second-order equation, 

 

 

Here unknown fields are now 𝑐 and 𝜇. The weak (variational) form of the problem reads:  

(𝑐, 𝜇) ∈ 𝑉×𝑉 such that 

 

 

 

 

C. Time discretization 

 

To deal any problem, the time derivative must be dealt with. To do this we apply the -method to the mixed weak form of equation: 

 

 

 

 
