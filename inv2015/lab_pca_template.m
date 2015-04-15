%
% lab_pca.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Principal Component Analysis of diets in Europe.
%

clc, clear, close all

% load the data, compute correlations, and make scatterplots

% get data (modified datafile from from protein.dat)
% X is the n x p predictor (or data) matrix
% n = number of observations (countries)
% p = number of variables (protein sources)
ddir = './data/';
X0 = load([ddir 'protein_matlab.dat']);      
[n,p] = size(X0);
xinds = [1:p];

% columns:
vars = {'RedMeat','WhiteMeat','Eggs','Milk','Fish','Cereals','Starch','Nuts','Fr&Veg'};
vlabs = {'RM','WM','EG','MK','FI','CL','ST','NT','FV'};

% rows:
country = {'Albania','Austria','Belgium','Bulgaria','Czechoslovakia',...
           'Denmark','E_Germany','Finland','France','Greece','Hungary',...
           'Ireland','Italy','Netherlands','Norway','Poland','Portugal',...
           'Romania','Spain','Sweden','Switzerland','UK','USSR','W_Germany','Yugoslavia'};
clabs =  {'ALB','AUS','BEL','BUL','CZK',...
           'DMK','EGR','FIN','FRA','GRE','HUN',...
           'IRE','ITY','NED','NOR','POL','POR',...
           'ROM','SPA','SWE','STZ','UK','USSR','WGR','YUG'};

whos      
    
% scatterplot matrix: uncentered, unstandardized
figure, gplotmatrix(X0,[],[],[],[],[],[],[],vars,vars), 

% START YOUR ANALYSIS HERE


