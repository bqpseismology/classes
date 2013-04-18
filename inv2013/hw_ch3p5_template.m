%
% hw_ch3p5.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Estimation of density variations using singular value decomposition.
%

clear
close all
clc

%==================================================================

% Exercise 3.5; TSVD problem

% Load the data (d)
load('/usr/local/matlab_toolboxes/aster/cd_5.2/Exercises/chap3/prob5/ifk.mat');
whos
m = length(d);  % number of data
n = m;          % number of model parameters

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

