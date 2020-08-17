%%
clear all, close all, clc, format long, format compact
%%
% *Exercise 2*
%
% (a)

% spectral radii
rhoBJ=[];
rhoBGS=[];
for N=[5,10,20,40,80]
    h = 1/N;
    A = (2/h^2)*diag(ones(N-1,1)) - (1/h^2)*diag(ones(N-2,1),1) - (1/h^2)*diag(ones(N-2,1),-1);
    b = sin(pi*h*(1:N-1)');
    D = diag(diag(A));
    E = D-tril(A);
    F = D-triu(A);
    BJ = D\(E+F);
    BGS = (D-E)\F;
    rhobj = max(abs(eig(BJ)));
    rhoBJ=[rhoBJ rhobj];
    rhobgs = max(abs(eig(BGS)));
    rhoBGS = [rhoBGS rhobgs];
end
rhoBJ'
rhoBGS'

% Both iterative methods should converge for these values of N as rho is <1 for each method Gauss-Seidel should converge faster as it has lower values of rho for each N. As N increases, rho tends to 1 therefore the methods should converge slower.

% Plotting 1-rho
hold off
Nvec=[5,10,20,40,80];
loglog(Nvec, 1-rhoBJ, 'b', 'LineWidth', 2)
hold on
loglog(Nvec, 1-rhoBGS, 'r', 'LineWidth', 2)
grid on
xlabel('N')
ylabel('1-rho')
legend('Jacobi', 'Gauss-Seidel', 'Location', 'SouthWest')

% eqn of the BJ line
mJ = (log(1-rhoBJ(1))-log(1-rhoBJ(5)))/(log(5)-log(80));
cJ = log(1 - rhoBJ(5)) - mJ*log(Nvec(5));
CJ = exp(cJ)
alphaJ = mJ

% eqn of BGS line
mGS = (log(1-rhoBGS(1))-log(1-rhoBGS(5)))/(log(5)-log(80));
cGS = log(1 - rhoBGS(5)) - mGS*log(Nvec(5));
CGS = exp(cGS)
alphaGS = mGS

% So the approximate relationship between the spectral radii and N is rho = 1 - C*N^alpha where C and alpha are defined above for each method.

% As N tends to infinity, 1-rho tends to 0 i.e rho tends to 1, so the number of iterations grows exponentially. Thus, when N is large, the methods will converge exponentially slower.
%%
% (b)

tol = 1e-10;
nmax = 1e5;

hold off
% Plot approximations
for N=[5,10,20,40,80]
    xvals=linspace(0, 1, N+1);
    x0 = zeros(N-1, 1);
    h = 1/N;
    A = (2/h^2)*diag(ones(N-1,1)) - (1/h^2)*diag(ones(N-2,1),1) - (1/h^2)*diag(ones(N-2,1),-1);
    b = sin(pi*h*(1:N-1)');
    [x, niter, relresiter, xiter] = itermeth(A, b, x0, nmax, tol, 'G');
    plot(xvals, [0;x;0])
    hold on
end
y=@(x)pi^(-2)*sin(pi*x);
xvals=linspace(0, 1, 81);
plot(xvals, y(xvals), 'black')

% As N increases, the solutions seem to converge.

%%
error_vect = [];
for N=[5,10,20,40,80]
    xvals=linspace(1/N, (N-1)/N, N-1);
    x0 = zeros(N-1, 1);
    h = 1/N;
    A = (2/h^2)*diag(ones(N-1,1)) - (1/h^2)*diag(ones(N-2,1),1) - (1/h^2)*diag(ones(N-2,1),-1);
    b = sin(pi*h*(1:N-1)');
    [x, niter, relresiter, xiter] = itermeth(A, b, x0, nmax, tol, 'G');
    diff=[];
    for i=1:N-1
        diff=[diff (x(i)-y(xvals(i)))];      
    end
    error_vect = [error_vect max(abs(diff))];
end
error_vect'

% Plotting N against the error
hold off
Nvec = [5,10,20,40,80];
hvect = 1./Nvec;
m = (log(error_vect(1))-log(error_vect(5)))/(log(hvect(1))-log(hvect(5)));
c = log(error_vect(5)) - m*log(hvect(5));
C = exp(c)
p = m
loglog(hvect, error_vect, 'blue')
% loglog(hvect, C*hvect.^m, 'red') (line check)

% My results support the theoretical estimate. We can see this because p = ~2.
%%
% (c)

tol = 1e-10;
nmax = 1e5;

% Error
y=@(x)(1/2)*x*(1-x);
eN = [];
for N=[5,10,20,40,80]
    xvals=linspace(1/N, (N-1)/N, N-1);
    x0 = zeros(N-1, 1);
    h = 1/N;
    A = (2/h^2)*diag(ones(N-1,1)) - (1/h^2)*diag(ones(N-2,1),1) - (1/h^2)*diag(ones(N-2,1),-1);
    b = ones(N-1,1);
    [x, niter, relresiter, xiter] = itermeth(A, b, x0, nmax, tol, 'G');
    diff=[];
    for i=1:N-1
        diff=[diff (x(i)-y(xvals(i)))];      
    end
    eN = [eN max(abs(diff))];
end

% Plotting the error
Nvec = [5,10,20,40,80];
hvect = 1./Nvec;
m = (log(eN(1))-log(eN(5)))/(log(hvect(1))-log(hvect(5)));
c = log(eN(5)) - m*log(hvect(5));
C = exp(c);
loglog(hvect, eN, 'blue')
% loglog(hvect, C*hvect.^m, 'red') (line check)

% It's surprising that the error should be so small. This has occured because C is proportional to y''''. But y'''' = 0 so we have C = 0. Now, as the error <= Ch^2, we should expect the error to be (close to) 0.
