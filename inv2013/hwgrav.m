%
% hwgrav.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Estimation of density from gravity using singular value decomposition.
%

clear, close all, clc

% directory for input data
dir0 = './';
ddir = [dir0 'data/'];
% directory for output figures
pdir = './';
iprint = 0;   % =1 to print figures to file

% get data and discretized integral
ww1 = 'gravity.dat';
ww2 = 'integral.dat';
load([ddir ww1]);
load([ddir ww2]);
g = gravity;
X = integral;

whos

% START YOUR WORK HERE


