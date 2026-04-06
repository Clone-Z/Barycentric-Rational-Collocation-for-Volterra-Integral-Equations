%% Example 3: Barycentric rational collocation for Volterra integral equation
%   This script solves the VIE:
%     u(x) - ∫_0^x (x-s)^{-0.4} s^{-0.2} u(s) ds = f(x),  0 ≤ x ≤ 1
%   with exact solution u(x) = x^{2/5} = x^{0.4}, and f(x) derived accordingly.
%   The method uses the 2nd-tapered BIEP rational interpolant (optimized for α=0.4)
%   and Gauss-Laguerre quadrature. The figure shows the root-exponential convergence
%   rate and the growth of the inverse norm of the discrete system.
%
%   References:
%     [33] Zhao & Xiang, IMA J. Numer. Anal., 2025
%     [29] Stahl, Bull. Amer. Math. Soc., 1993
%
%   Author: Based on code from Zhong, Zhao & Xiang (2026)
%   Modified: with comments for GitHub

clear; clc;
tic
%% Load precomputed barycentric weights and nodes (optimized for α=0.4)
%   The file 'NW04.mat' contains cell arrays NW{n,1} (nodes) and NW{n,2} (weights)
%   for the rational interpolation with tapered exponentially clustered poles.
%   Here n corresponds to the number of interpolation points (10,20,...,240).
load NW04.mat

%% Define a fine evaluation grid for error calculation
xx = [exp(chebpts(1e4, [-300, -10])); chebpts(1e4, [exp(-10), 1])];  % 20,000 points

%% Parameter settings
nn = 10:10:240;                     % numbers of interpolation points
k = 1;                              % counter for error storage
invNorm = zeros(size(nn));          % infinity norm of the inverse matrix (I - A)^{-1}
err = zeros(size(nn));              % uniform-norm error on [0,1]

%% Problem data
%   Exact solution and right-hand side
%   Kernel has diagonal and boundary singularities: (x-s)^{-0.4} s^{-0.2}
f  = @(t) t.^(0.4) - t.^(0.8) * gamma(3/5) * gamma(6/5) / gamma(9/5);
ut = @(t) t.^(2/5);                 % exact solution: u(x) = x^{0.4}
Ut = ut(xx);                        % exact solution on fine grid

%   Weak singularity parameters
alp  = -0.4;                        % exponent of (x-s) in kernel
beta = -0.2;                        % exponent of s in kernel

%   Function used in the transformed integral (Equation (20))
%   q(s,t) = t^{1+alp+beta} * ((1-exp(-s))/s)^{alp} * exp(-beta*s)
q = @(s, t) t.^(1 + alp + beta) .* ((1 - exp(-s)) ./ s).^alp .* exp(-beta * s);

%% Main loop over increasing n
for n = nn
    % 1) Gauss-Laguerre nodes and weights for the integral transformation
    %    Use lagpts (or glags) to obtain nodes and weights for weight y^{-λ} e^{-y}
    [X, W] = lagpts(n + 1, alp);    % alp = -0.4, so weight y^{0.4} e^{-y}
    %    Alternative: [X, W] = glags(n+5, alp); (if glags is available)

    % 2) Retrieve BIEP nodes and weights for rational interpolation
    xi = NW{n/10, 1};               % interpolation nodes (size n+1)
    wi = NW{n/10, 2};               % barycentric weights for these nodes

    % 3) Assemble the discretized system matrix L (size (n+1)x(n+1))
    %    L(i,j) = x_i^{1-λ-θ} * ∫_0^∞ k(x_i, s) l_j(x_i e^{-s/(1-θ)}) s^{-λ} e^{-s} ds
    %    Here λ = -alp, θ = -beta, so 1-λ-θ = 1+alp+beta.
    E = eye(n+1);
    L1 = zeros(size(xx));       % accumulate over fine grid points
    for i = 1:length(xi)
        L1 = L1 + abs(bary(xx, E(:,i), xi, wi));
    end
    LC1(k) = max(L1);
    L = zeros(n+1);
    for i = 1:(n+1)
        t = xi(i);
        for j = 1:(n+1)
            S1 = t * exp(-X);       % s = x_i * exp(-y) in transformation
            int1 = W * (q(X, t) .* bary(S1, E(:,j), xi, wi));
            L(i,j) = int1;
        end
    end
      LC1(k) = max(L1);
    % 4) Build the linear system (I - L) u_n = f
    %    The original equation is u - integral = f.
    N = E - L;
    F = f(xi);                      % right-hand side at interpolation nodes

    %    Compute the infinity norm of the inverse of the system matrix
    invNorm(k) = norm(inv(N), 'inf');

    %    Solve for the approximate solution at nodes
    Ui = N \ F;

    % 5) Evaluate the rational interpolant on the fine grid and compute uniform error
    U = bary(xx, Ui, xi, wi);
    err(k) = norm(U - Ut, 'inf');

    k = k + 1;
end

%% Plot results: convergence error (left) and inverse norm (right)
figure;
tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact');

% Left subplot: uniform-norm error vs sqrt(n)
nexttile;
semilogy(sqrt(nn), err, 'o-', 'LineWidth', 1.2); hold on;
% Theoretical optimal rate (Stahl) for α=0.4
semilogy(3:10, 1 * exp(-2 * pi * sqrt(0.4) * (3:10)), '--k', 'LineWidth', 1.2);
grid on;
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times');
xlabel('$\sqrt{n}$', 'Interpreter', 'latex');
ylabel('uniform-norm error');
legend('Errors on $[0,1]$', '$100\exp(-2\pi\sqrt{0.4 n})$', ...
       'Interpreter', 'latex');
axis([2, 16, 1e-16, 1e-2]);

% Right subplot: growth of the inverse norm ||(I-A)^{-1}||_∞
% nexttile;
% semilogy(nn, invNorm, 's-', 'LineWidth', 1.2);
% grid on;
% set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times');
% xlabel('$n$', 'Interpreter', 'latex');
% %ylabel('$\|(I-\hat{A})^{-1}\|_{\infty}$', 'Interpreter', 'latex');
nexttile
semilogx(nn, invNorm, '*-', 'LineWidth', 1.2); hold on
semilogx(nn, LC1,  'o-', 'LineWidth', 1.2); grid on
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
legend('$\|(I-\hat{A})^{-1}\|_{\infty}$', 'Lebesgue constant', 'Interpreter', 'latex')
xlabel('${n}$', 'Interpreter', 'latex')
axis([0, 300, 0, 35]);

title('Stability measure');
%axis([2, 240, 20, 30])
% Optional: Also plot pointwise error for the largest n (n=240) if desired
% figure;
% semilogy(xx, abs(U - Ut), '-', 'LineWidth', 1.2);
% grid on;
% xlabel('$x$', 'Interpreter', 'latex');
% ylabel('absolute error');
% title('Pointwise error for n=240');
toc
