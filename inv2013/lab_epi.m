%
% lab_epi.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Example script for generating samples of a function using the rejection
% method (Tarantola 2005, Section 2.3.2).
%
% See notes for lab exercise lab_epi.pdf
%
% This script is a guide to the epicenter homework problem presented as
% Problem 7-1 in Tarantola (2005).
%
% calls plot_histo.m
%

clear, close all, clc

% limits for x to consider
% OUTSIDE THESE LIMITS WE ASSUME THAT f(x) = 0
xmin = -15; xmax = 12;

% define in-line function f(x)
x1 = -2; A1 = 2; sig1 = 2;
x2 =  4; A2 = 1; sig2 = 0.5;
% note that parameters x0 and A must be previously defined
f = @(x) ( A1*exp(-(x-x1).^2/(2*sig1^2)) );
%f = @(x) ( A1*exp(-(x-x1).^2/(2*sig1^2)) + A2*exp(-(x-x2).^2/(2*sig2^2)) );

% KEY TECHNICAL POINT: f is a function, not a numerical array
% (note that x is not a stored array)
whos

% generate samples
NTRY = 1e5;
xtry = xmin + (xmax-xmin)*rand(NTRY,1);

% sample the function
fmax = max([A1 A2]);
ftry = f(xtry) / fmax;          % SET A: values between 0 and 1
chance = rand(NTRY,1);          % SET B: values between 0 and 1

% KEY COMMAND: compare pairs of test samples in sets A and B,
%              then accept or reject the test sample
xkeep = xtry(ftry > chance);

% plot
xcurve = linspace(xmin,xmax,1000);
figure; nr=3; nc=2;
edges1 = [xmin:0.2:xmax];
edges2 = [0:0.05:1];
subplot(nr,nc,1); plot(xcurve,f(xcurve)/fmax);
xlabel('x'); ylabel('f(x)'); title('(a)'); axis([xmin xmax 0 1.2]);
subplot(nr,nc,2); plot_histo(xtry,edges1);
xlabel('xtry'); title('(b)'); 
subplot(nr,nc,3); plot_histo(ftry,edges2);
xlabel('f(xtry)'); title('(c)'); 
subplot(nr,nc,4); plot_histo(chance,edges2);
xlabel('chance'); title('(d)'); 
subplot(nr,nc,5); plot_histo(xkeep,edges1);
xlabel('xkeep');  title('(e)'); xlim([xmin xmax]);

%==========================================================================
