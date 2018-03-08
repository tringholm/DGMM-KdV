%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF USING THIS CODE FOR RESEARCH PURPOSES, PLEASE CITE OUR ARTICLE     %
% Eidnes, S., Owren, B. & Ringholm, T. Adv Comput Math (2017).          %
% https://doi.org/10.1007/s10444-017-9562-8                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial solution: u(x,0) = c/2*sech^2(sqrt(c)/2*x)
% Periodic boundary conditions: u(L,t) = u(-L,t)
% Preserving H = \int(.5*u'^2 - u^3)dx
% Using Legendre-Gauss quadrature to evaluate integrals, employing the code lgwt.m written by Greg von Winckel

clear
close all

M = 100;        % Number of spatial discretization points (-1)
k = .1;         % Time step
T = 10;         % End time
L = 100;        % Interval [-L,L]
moving = true;  % true if using moving mesh
doplot = true;  % true if producing plots while running
method = 1;     % 1 for AVF DG method, 2 for implicit midpoint method
order = 3;      % order of the basis functions
monitor = 1;    % 1 for arc-length monitor function, 2 for curvature

% Analytic solution
c = 6; 
u_analytic = @(x,t) .5*c*sech(0.5*sqrt(c)*min(abs(x-c*t),abs(2*L+x-c*t))).^2;

% Run function with plotting of analytic solutions together with numerical
% solution
[u,p,H] = KdV_AVF(M,k,T,L,method,order,moving,monitor,doplot,u_analytic);