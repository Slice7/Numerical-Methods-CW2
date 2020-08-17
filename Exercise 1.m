%%
clear all, close all, clc, format long, format compact
%%  
% *Exercise 1*
% 
% (a)

% To show A is SDD we need to solve 2|epsilon^2| + 2|epsilon| < 1. Consider f(epsilon) = 2|epsilon^2| + 2|epsilon| - 1 = 0. Solving this gives epsilon = +/-(-1+sqrt(3))/2, so for f < 0, we have -(-1+sqrt(3))/2 < epsilon < (-1+sqrt(3))/2. Back to the original question, as epsilon is contained in [0, 1], we have 0 < epsilon < (-1+sqrt(3))/2.
%%
% (b)

n = 5;
epsi = 0.3;
tol=1e-10;
nmax=1e3;
x0=[0,0,0,0,0]';
[A, b] = matrix(n, epsi);
% Jacobi method
[x, niter, relresiter, xiter] = itermeth(A, b, x0, nmax, tol, 'J');
% Gauss-Seidel method
[x, niter, relresiter, xiter] = itermeth(A, b, x0, nmax, tol, 'G');

% We see that the Jacobi method converged in 50 iterations and the Gauss-Seidel method converged in 14 iterations.
%%
% (c)

specBJ = [];
specBGS = [];
for i=0:100
    [A, b] = matrix(n, i/100);
    D = diag(diag(A));
    E = D-tril(A);
    F = D-triu(A);
    BJ = inv(D)*(E+F);
    specBJ = [specBJ max(abs(eig(BJ)))];
    BGS = inv(D-E)*F;
    specBGS = [specBGS max(abs(eig(BGS)))];
end
% plot
x = linspace(0, 1, 101);
plot(x, specBJ, 'b', 'LineWidth', 2)
hold
plot(x, specBGS, 'r', 'LineWidth', 2)
grid on
xlabel('epsi')
ylabel('rho')
legend('Jacobi', 'Gauss-Seidel', 'Location', 'NorthWest')

% The plots tell us that as epsilon grows both methods converge slower (as the spectral radii get bigger.) We'd expect Jacobi to converge for epsilon < ~0.44 and Gauss-Seidel to converge for epsilon < ~0.7. This is a wider range than in part (a). This is expected because A being SDD is only a sufficient condition.
%%
% (d)

% For n=5 and epsilon=0.5, Gauss-Seidel would be recommended because its iteration matrices have lower spectral radii compared to the Jacobi method, therefore we will have faster convergence.
