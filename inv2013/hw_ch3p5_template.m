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

m = length(d);

% set up vector of the y collocation points
dx = 0.05;
xmin = dx/2;
xmax = 1 - dx/2;
y = (xmin:dx:xmax)';

% set up vectors of the x collocation points


