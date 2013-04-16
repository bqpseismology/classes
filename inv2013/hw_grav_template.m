%
% hw_grav.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Estimate density variations using singular value decomposition.
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
d = gravity;
n = length(d);
m = n;

% construct discretization vector for data
xmin = 0;
xmax = 1;
xvec = collocate(xmin,xmax,m);

% construct discretization vector for model
ymin = 0;
ymax = 1;
yvec = collocate(ymin,ymax,n);

whos xvec yvec d

% example of using plotconst_mod.m
figure; hold on;
plotconst_mod(rand(m,1),xmin,xmax,{'k','linewidth',2});
plotconst_mod(rand(m,1),xmin,xmax,{'r--','linewidth',2});

% construct design matrix G

