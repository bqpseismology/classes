%
% newton.m
% Carl Tape, GEOS 694, Applied Seismology, Feb 2012
%
% This is a template for implementing an iterative netwon algorithm.
% See Tarantola (2005), Eq. 6.291.
%
% LAB EXERCISE:
%   1. Using 'fminbnd' compute the (numerical) minimum of F(m).
%       Plot the point (mmin,F(mmin)) on the curve F(m).
%   2. Implement the Newton algorithm and show that you reach the minimum.
%       Plot the point (m,F(m)) for each new value of m.
%       Try severl different starting values of m to make sure it works.
%   3. Define a sensible stopping criterion.
%   4. Repeat the experiment for a quadratic function P = [2 -1 2].
% 
%

close all
clear
clc
format long

% define a misfit function and its derivatives
P = [2 -1 2 -1 2];
%F = @(m) ( polyval(P,m) );
F   = @(m) ( P(1)*m.^4 + P(2)*m.^3 + P(3)*m.^2 + P(4)*m + P(5) );
dF1 = @(m) ( 4*P(1)*m.^3 + 3*P(2)*m.^2 + 2*P(3)*m + P(4) );
dF2 = @(m) ( 12*P(1)*m.^2 + 6*P(2)*m + 2*P(3) );

% specify bounds for choosing initial model (and for plotting)
mA = -2;
mB = 2;
m = linspace(mA,mB,100);

figure; hold on;
plot(m,F(m));
xlabel('model m'); ylabel('misfit function, F(m)');

% initial model
%m0 = -1.5;
m0 = mA + (mB-mA)*rand;

%==========================================================================
