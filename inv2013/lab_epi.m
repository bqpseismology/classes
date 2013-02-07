%
% lab_epi.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Example script for using an in-line function, then generating samples as
% if the function were a probability distribution.
%
% The sampling method here is called the Rejection Method (Tarantola 2005, Section 2.3.2),
% which is attributed to a paper by John von Neumann (1951).
% There is a function in the Matlab Statistics Toolbox (accrejrnd):
%    http://www.mathworks.com/help/stats/common-generation-methods.html#br5k9hi-4
%
% This script is a guide to the epicenter homework problem presented as
% Problem 7-1 in Tarantola (2005).
%
% calls plot_histo.m
%

clear, close all, clc

% limits for x to consider
xmin = -10; xmax = 10;

% in-line function f(x)
x0 = -2;
A = 2;
% note that parameters x0 and A must be previously defined
f = @(x) ( A*exp(-(x-x0).^2) );

% KEY TECHNICAL POINT: f is a function, not a numerical array
% (note that x is not a stored array)
whos

% generate samples
NTRY = 1e5;
xtry = xmin + (xmax-xmin)*rand(NTRY,1);

% sample the function
fmax = A;
ftry = f(xtry) / fmax;      % value between 0 and 1
chance = rand(NTRY,1);
xkeep = xtry(chance < ftry);

% plot
xcurve = linspace(xmin,xmax,100);
figure; nr=3; nc=2;
edges1 = [xmin:0.2:xmax];
edges2 = [0:0.05:1];
subplot(nr,nc,1); plot(xcurve,f(xcurve)/fmax); ylim([0 1]);
subplot(nr,nc,2); plot_histo(xtry,edges1); xlabel('xtry');
subplot(nr,nc,3); plot_histo(ftry,edges2); xlabel('f(xtry)');
subplot(nr,nc,4); plot_histo(chance,edges2); xlabel('chance');
subplot(nr,nc,5); plot_histo(xkeep,edges1); xlabel('xkeep');

% START EXERCISE HERE


%==========================================================================
