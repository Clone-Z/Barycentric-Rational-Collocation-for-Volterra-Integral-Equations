%% Example 1: Barycentric rational collocation for Volterra integral equation
%   This script solves the VIE:
%     u(x) + (1/sqrt(pi)) * ∫_0^x (x-s)^{-1/2} s^{-1/2} u(s) ds = f(x)
%   with exact solution u(x) = x^{1/pi}, and f(x) defined accordingly.
%   The method uses the 2nd-tapered BIEP rational interpolant and Gauss-Laguerre
%   quadrature for the integrals. The convergence is compared with Floater-Hormann
%   rational interpolation.
%
%   References:
%     [33] Zhao & Xiang, IMA J. Numer. Anal., 2025
%     [29] Stahl, Bull. Amer. Math. Soc., 1993
%
%   Author: Based on code from Zhong, Zhao & Xiang (2026)
%   Modified: with comments for GitHub

clear; clc;

%% Load precomputed barycentric weights and nodes for 2nd-tapered BIEP
%   The file 'NWpi.mat' contains cell arrays NW{n,1} (nodes) and NW{n,2} (weights)
%   for the rational interpolation with tapered exponentially clustered poles.
%   Here n corresponds to the number of interpolation points (10,20,...,240).
load NWpi.mat

%% Define a fine evaluation grid for error calculation
%   Use Chebyshev points on two intervals to cover [0,1] accurately.
xx = [exp(chebpts(1e4, [-300, -10])); chebpts(1e4, [exp(-10), 1])];  % 20,000 points

%% Parameter settings
nn = 10:10:240;             % numbers of interpolation points
k = 1;                      % counter for error storage
err1 = zeros(length(xx), 3);  % absolute errors for selected n (30,100,240)
err2 = zeros(size(nn));        % uniform-norm errors for all n
invNorm = zeros(size(nn));     % norm of inverse of the discrete system matrix (I - A)^{-1}
LC1  = zeros(size(nn));        % Lebesgue constant
mm = [30 100 240];              % indices for which we store pointwise errors

%% Problem data
%   Exact solution and right-hand side
%   alpha = 1/pi, lambda = theta = 0.5, so the solution has a weak singularity.
f  = @(t) (1 + gamma(1/2 + 1/pi) / gamma(1 + 1/pi)) * t.^(1/pi);
ut = @(t) t.^(1/pi);
Ut = ut(xx);                     % exact solution on fine grid

%   Weak singularity parameters (the kernel is (x-s)^{-1/2} s^{-1/2})
alp = -0.5;    % exponent of (x-s) in kernel
beta = -0.5;   % exponent of s in kernel

%   Function used in the transformed integral (Equation (20))
%   q(s,t) = t^{1+alp+beta} * ((1-exp(-s))/s)^{alp} * exp(-beta*s)
q = @(s, t) t.^(1 + alp + beta) .* ((1 - exp(-s)) ./ s).^alp .* exp(-beta * s);

%% Main loop over increasing n
for n = nn
    % 1) Gauss-Laguerre nodes and weights for the integral transformation
    %    Use lagpts (or glags) to obtain nodes and weights for weight y^{-λ} e^{-y}
    %    Here we take N = n+20 points (more than n to ensure high quadrature accuracy)
    [X, W] = lagpts(n + 1, alp);   % alp = -0.5, so weight y^{0.5} e^{-y}
    

    % 2) Retrieve BIEP nodes and weights for rational interpolation
    xi = NW{n/10, 1};   % interpolation nodes (size n+1)
    wi = NW{n/10, 2};   % barycentric weights for these nodes

    % 3) Compute Lebesgue constant for the interpolation operator
    %    Λ_n = max_{x∈[0,1]} Σ_j |l_j(x)|
    E = eye(n+1);               % identity matrix, each column is a unit vector
    L1 = zeros(size(xx));       % accumulate over fine grid points
    for i = 1:length(xi)
        L1 = L1 + abs(bary(xx, E(:,i), xi, wi));
    end
    LC1(k) = max(L1);

    % 4) Assemble the discretized system matrix L (size (n+1)x(n+1))
    %    L(i,j) = x_i^{1-λ-θ} * ∫_0^∞ k(x_i, s) l_j(x_i e^{-s/(1-θ)}) s^{-λ} e^{-s} ds
    %    Here λ = -alp, θ = -beta, so 1-λ-θ = 1+alp+beta = 1 -0.5 -0.5 = 0.
    L = zeros(n+1);
    for i = 1:(n+1)
        t = xi(i);              % current interpolation point
        for j = 1:(n+1)
            % evaluate integrand at Gauss-Laguerre nodes
            S1 = t * exp(-X);   % s = x_i * exp(-y) in transformation
            % compute the barycentric rational basis l_j at points S1
            int1 = W * (q(X, t) .* bary(S1, E(:,j), xi, wi));
            L(i,j) = int1;
        end
    end

    % 5) Build the linear system (I + (1/sqrt(pi)) L) u_n = f
    %    Since the original equation is u + (1/sqrt(pi)) * integral = f.
    N = E + (1/sqrt(pi)) * L;
    F = f(xi);                  % right-hand side at interpolation nodes

    %    Compute the infinity norm of the inverse of the system matrix
    %    This measures the stability of the discrete system.
    invNorm(k) = norm(inv(N), 'inf');

    %    Solve for the approximate solution at nodes
    Ui = N \ F;

    % 6) Evaluate the rational interpolant on the fine grid
    U = bary(xx, Ui, xi, wi);

    % 7) Store errors
    if any(mm == n)
        err1(:, mm == n) = abs(U - Ut);
    end
    err2(k) = norm(U - Ut, 'inf');
    k = k + 1;
end

%% Plot results
%   Left subplot: uniform-norm error vs sqrt(n)
%   Right subplot: pointwise absolute error for selected n

tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact')

nexttile
semilogy(sqrt(nn), err2, 'o-', 'LineWidth', 1.2); hold on
% theoretical optimal rate (Stahl)
semilogy(3:10, 10 * exp(-2 * pi * sqrt(1/pi) * (3:10)), '--k', 'LineWidth', 1.2);
% mark points for n = 30,100,240
semilogy(sqrt(30),  err2(3),  '.k', 'MarkerSize', 12);
semilogy(sqrt(100), err2(10), '.k', 'MarkerSize', 12);
semilogy(sqrt(240), err2(24), '.k', 'MarkerSize', 12);
% load precomputed Floater-Hormann results (if available)
load Figure_8FH5.mat   % contains err2_FH for n=10:10:240
semilogy(sqrt(nn), err2(1:24), '*-k', 'LineWidth', 1); hold on
% algebraic convergence of Floater-Hormann
semilogy(sqrt(nn), nn.^(-1/pi) / 10, ':k', 'LineWidth', 1.2);
grid on
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
xlabel('$\sqrt{n}$', 'Interpreter', 'latex')
ylabel('uniform-norm error')
legend('$\alpha=1/\pi$', '$10\exp(-2\pi\sqrt{n/\pi})$', '', '', '', ...
       'F-H ($d=5$)', '$n^{-1/\pi}/10$', 'Interpreter', 'latex')
axis([2, 16, 1e-16, 1e0])

nexttile
loglog(xx, err1(:,1), '-',  'LineWidth', 1.2); hold on
loglog(xx, err1(:,2), '-',  'LineWidth', 1.2);
loglog(xx, abs(U - Ut), '-', 'LineWidth', 1.2);  % n=240
grid on
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
legend('$n=30$', '$n=100$', '$n=240$', 'Interpreter', 'latex')
xlabel('$x$', 'Interpreter', 'latex')
ylabel('absolute error')
axis([1e-100, 1e0, 1e-40, 1e0])

%% Additional figure: norm of inverse and Lebesgue constant
figure(2)
tiledlayout(1, 2, 'TileSpacing', 'loose', 'Padding', 'compact')

nexttile
plot(nn, invNorm, '*-', 'LineWidth', 1.2); hold on
plot(nn, LC1,  'o-', 'LineWidth', 1.2); grid on
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
legend('$\|(I-\hat{A})^{-1}\|_{\infty}$', 'Lebesgue constant', 'Interpreter', 'latex')
xlabel('${n}$', 'Interpreter', 'latex')

nexttile
semilogx(nn, invNorm, '*-', 'LineWidth', 1.2); hold on
semilogx(nn, LC1,  'o-', 'LineWidth', 1.2); grid on
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
legend('$\|(I-\hat{A})^{-1}\|_{\infty}$', 'Lebesgue constant', 'Interpreter', 'latex')
xlabel('${n}$', 'Interpreter', 'latex')