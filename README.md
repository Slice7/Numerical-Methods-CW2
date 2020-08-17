# Exercise 1
Consider the *n* x *n* band matrix *A*<sub>*ϵ*</sub>, whose values are 1 on the main diagonal, ϵ on the first diagonals, and ϵ<sup>2</sup> on the second diagonals, ϵ ∈ [0, 1]. For example, if *n* = 5, we'd have the following matrix:
|        1        |        ϵ        |  ϵ<sup>2</sup>  |        0        |        0        |
|:---------------:|:---------------:|:---------------:|:---------------:|:---------------:|
|      **ϵ**      |      **1**      |      **ϵ**      |**ϵ<sup>2</sup>**|      **0**      |
|**ϵ<sup>2</sup>**|      **ϵ**      |      **1**      |      **ϵ**      |**ϵ<sup>2</sup>**|
|      **0**      |**ϵ<sup>2</sup>**|      **ϵ**      |      **1**      |      **ϵ**      |
|      **0**      |      **0**      |**ϵ<sup>2</sup>**|      **ϵ**      |      **1**      |

We'll be looking at the linear system *A*<sub>*ϵ*</sub>**x** = **b**<sub>*ϵ*</sub>, where **b**<sub>*ϵ*</sub> = *A*<sub>*ϵ*</sub>(1, 1, **...** , 1)<sup>*T*</sup>.

**(a)** First, we find the values of *ϵ* for which *A*<sub>*ϵ*</sub> is strictly diagonally dominant.

**(b)** Now we'll solve the system for *n* = 5 and *ϵ* = 0.3, using both the Jacobi and Gauss-Seidel methods, noting the convergence rate for each method.

**(c)** We then compute the spectral radius *ρ* of the iteration matrices of each method for *ϵ* = 0, 0.01, 0.02, **...** , 1, and plot the results onto a graph. This should tell us how the two methods perform for different values of *ϵ*.

**(d)** Lastly, we conclude which method would be recommended for *n* = 5 and *ϵ* = 0.5, based on the graph from **(c)**.

# Exercise 2
Consider the boundary value problem

- *d*<sup>2</sup>*y*/*dx*<sup>2</sup> = -sin(π*x*), *x* ∈ [0, 1],  
- *y*(0) = *y*(1) = 0.

The exact solution is *y*(*x*) = π<sup>-2</sup>sin(π*x*).

The solution can be approximated on a set of *N*+1 equispaced nodes {*x*<sub>0</sub>, *x*<sub>1</sub>, **...** ,*x*<sub>*N*</sub>}, where *x*<sub>0</sub> = 0 and *x*<sub>*N*</sub> = 1, *h* = 1/*N* is the step size, and *y*(*x*<sub>*n*</sub>) is approximated by *u*<sub>*n*</sub>.

We obtain the approximation **u** = (*u*<sub>1</sub>, **...** , *u*<sub>*N*-1</sub>)<sup>*T*</sup>, which satisfies the equation *A***u** = **b**.

**(a)** We will again compute the spectral radius of the iteration matrices of the Jacobi and Gauss-Seidel methods for the system above, for *N* ∈ {5, 10, 20, 40, 80}, and analyse the results with respect to the performance of each method.

**(b)** We then solve the system using the Gauss-Seidel method and plot the 5 approximations on the same axes, along with the exact solution for comparison.  
On a separate graph, we'll plot the maximum error of each approximation against the corresponding step size, and examine the results.

**(c)** Now we consider the same boundary value problem, but this time with *d*<sup>2</sup>*y*/*dx*<sup>2</sup> = -1.  
Once again, we'll plot the maximum error of 5 approximations against their step sizes (using the same values of *N*), and comment on our findings.

# Concepts
Here are some brief explanations of the concepts used in these exercises.

### Strictly diagonally dominant
A matrix is strictly diagonally dominant if, for any given row in the matrix, the diagonal value has larger magnitude than the sum of the magnitudes of all other values in that row.

If a matrix is SDD, then both the Jacobi and Gauss-Seidel methods converge to the solution.

### Spectral radius
The spectral radius of a matrix *B* is simply its largest eigenvalue by absolute value. We denote this *ρ*(*B*).  
If *ρ*(*B*) < 1, the error **x** - **x**<sub>n</sub> → 0 as n → ∞. The converse is also true, but not relevant to these exercises.

### Jacobi method
This method decomposes an *n* x *n* matrix *A* into three matrices *D* - (*L* + *U*) = *A*,  
where *D* is the diagonal part of *A*, *L* is the strictly lower triangular part of -*A*, and *U* is the strictly upper triangular part of -*A*.  
This gives an iteration matrix *B<sub>J</sub>* = *D*<sup>-1</sup>(*L* + *U*)

We then solve the system using the iteration **x**<sub>*n*+1</sub> = *B*<sub>*J*</sub>**x**<sub>*n*</sub> + **c**, where **c** a vector derived from **b**.

### Gauss-Seidel method
This method uses the same decomposition as the Jacobi method, but uses a different iteration matrix, namely  
*B<sub>GS</sub>* = (*D* - *L*)<sup>-1</sup>*U*
