%% Example 2: Barycentric rational collocation for Volterra integral equation
%   This script solves the VIE:
%     u(x) + ∫_0^x (x-s)^{-0.2} u(s) ds = f(x),  0 ≤ x ≤ 1
%   with exact solution u(x) = x^{1/5} (alpha = 0.2), and f(x) derived accordingly.
%   The method uses the 2nd-tapered BIEP rational interpolant and Gauss-Laguerre
%   quadrature for the integrals. The convergence is compared with the optimal
%   root-exponential rate.
%
%   References:
%     [33] Zhao & Xiang, IMA J. Numer. Anal., 2025
%     [29] Stahl, Bull. Amer. Math. Soc., 1993
%
%   Author: Based on code from Zhong, Zhao & Xiang (2026)
%   Modified: with comments for GitHub

clear; clc;

%% Load precomputed barycentric weights and nodes for 2nd-tapered BIEP
%   The file 'NW02.mat' contains cell arrays NW{n,1} (nodes) and NW{n,2} (weights)
%   for the rational interpolation with tapered exponentially clustered poles.
%   Here n corresponds to the number of interpolation points (10,20,...,240).
%   The nodes and weights are optimized for approximating functions with branch
%   singularity x^{0.2} on [0,1].
load NW02.mat

%% Define a fine evaluation grid for error calculation
%   Use Chebyshev points on two intervals to cover [0,1] accurately.
xx = [exp(chebpts(1e4, [-300, -10])); chebpts(1e4, [exp(-10), 1])];  % 20,000 points

%% Parameter settings
nn = 10:10:240;                     % numbers of interpolation points
k = 1;                              % counter for error storage
err1 = zeros(size(nn));             % errors at collocation points
err2 = zeros(size(nn));             % uniform-norm errors on [0,1]

%% Problem data
%   Exact solution and right-hand side
%   The kernel has a weak diagonal singularity: (x-s)^{-0.2} (theta = 0)
%   The solution u(x) = x^{1/5} = x^{0.2}
f  = @(t) (pi/5) * t .* csc(pi/5) + t.^(1/5);
ut = @(t) t.^(1/5);
Ut = ut(xx);                        % exact solution on fine grid

%   Weak singularity parameter (lambda = 0.2)
alp = -0.2;                         % exponent of (x-s) in kernel

%   Function used in the transformed integral (Equation (20))
%   For theta = 0, the transformation simplifies: q(s,t) = t^{1+alp} * ((1-exp(-s))/s)^{alp}
q = @(s, t) t.^(1 + alp) .* ((1 - exp(-s)) ./ s).^alp;

%% Main loop over increasing n
for n = nn
    % 1) Gauss-Laguerre nodes and weights for the integral transformation
    %    Use lagpts (or glags) to obtain nodes and weights for weight y^{-λ} e^{-y}
    %    Here we take N = n+1 points (more than n to ensure high quadrature accuracy)
    [X, W] = lagpts(n + 11, alp);    % alp = -0.2, so weight y^{0.2} e^{-y}
    %    Alternative: [X, W] = glags(n+1, alp); (if glags is available)

    % 2) Retrieve BIEP nodes and weights for rational interpolation
    xi = NW{n/10, 1};               % interpolation nodes (size n+1)
    wi = NW{n/10, 2};               % barycentric weights for these nodes

    % 3) Assemble the discretized system matrix L (size (n+1)x(n+1))
    %    L(i,j) = x_i^{1-λ} * ∫_0^∞ k(x_i, s) l_j(x_i e^{-s/(1-θ)}) s^{-λ} e^{-s} ds
    %    Here θ = 0, so the factor becomes x_i^{1-λ}, and the integration variable
    %    transformation uses s = x_i * exp(-y) directly.
    E = eye(n+1);                   % identity matrix
    L = zeros(n+1);
    for i = 1:(n+1)
        t = xi(i);                  % current interpolation point
        for j = 1:(n+1)
            S1 = t * exp(-X);       % s = x_i * exp(-y) in transformation
            % compute the barycentric rational basis l_j at points S1
            int1 = W * (q(X, t) .* bary(S1, E(:,j), xi, wi));
            L(i,j) = int1;
        end
    end

    % 4) Build the linear system (I + L) u_n = f
    %    Since the original equation is u + integral = f.
    N = E + L;
    F = f(xi);                      % right-hand side at interpolation nodes

    %    Solve for the approximate solution at nodes
    Ui = N \ F;

    % 5) Compute errors
    %    err1: error at collocation points (max norm)
    err1(k) = norm(Ui - ut(xi), 'inf');

    %    Evaluate the rational interpolant on the fine grid
    U = bary(xx, Ui, xi, wi);

    %    err2: uniform error on the fine grid
    err2(k) = norm(U - Ut, 'inf');

    k = k + 1;
end

%% Plot results
%   Left/upper subplot: uniform-norm errors vs sqrt(n)
%   Right/lower subplot: pointwise absolute error for n=240

figure(1)
semilogy(sqrt(nn), err1, '*-', 'LineWidth', 1.2); hold on
semilogy(sqrt(nn), err2, 'o-', 'LineWidth', 1.2);
% theoretical optimal rate (Stahl) for alpha = 0.2
semilogy(3:14, 100 * exp(-2 * pi * sqrt(0.2) * (3:14)), '--k', 'LineWidth', 1.2);
grid on
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
xlabel('$\sqrt{n}$', 'Interpreter', 'latex')
ylabel('uniform-norm error')
legend('Errors at collocation points', 'Errors on $[0,1]$', ...
       '$100\exp(-2\pi\sqrt{0.2 n})$', 'Interpreter', 'latex')
axis([2, 16, 1e-16, 1e-2])

% figure(2)
% subplot(2,1,1)
% semilogy(xx, abs(U - Ut), '-', 'LineWidth', 1.2); grid on
% set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
% xlabel('$x$', 'Interpreter', 'latex')
% ylabel('absolute error')
% axis([0, 1, 1e-40, 1e0])
% 
% subplot(2,1,2)
% loglog(xx, abs(U - Ut), '-', 'LineWidth', 1.2); grid on
% set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times')
% xlabel('$x$', 'Interpreter', 'latex')
% ylabel('absolute error')
% axis([1e-120, 1e0, 1e-80, 1e0])