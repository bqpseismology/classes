%
% hw_grav.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Estimation of density variations using singular value decomposition.
%

clear
close all
clc

% directory for input data
dir0 = './';
ddir = [dir0 'data/'];

%==========================================================================

% get data and discretized integral
load([ddir 'gravity.dat']);
load([ddir 'integral.dat']);
g = gravity;
X = integral;
whos g X

[m,n] = size(X);


