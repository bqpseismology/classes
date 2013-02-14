%
% hwcov.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Template for HW on covariance and Gaussian random fields.
%
% calls covC.m
%

clear, close all, clc

ddir = './data/';

% load data
% m_samples : M x P
%         x : M x 1
% M is the number of model parameters describing a single sample (of a covariance matrix)
% P is the number of samples
load([ddir 'covhwdata']);
whos
[M,P] = size(m_samples);
xmin = min(x);
xmax = max(x);
ymin = min(m_samples(:));
ymax = max(m_samples(:));
ax0 = [xmin-1 xmax+1 ymin ymax];    % axes for plotting samples

% compute grid of distances among all points
% (D is needed for calculating the covariance matrix)
[X1,X2] = meshgrid(x,x);
D = abs(X1-X2);

% example of plotting with imagesc
figure; imagesc(D); xlabel('x index'); ylabel('x index');
title('distance between points'); colorbar

% START HERE




break

% CODE FOR COMPUTING NORMS OF YOUR ESTIMATED SAMPLES
% THIS ASSUMES YOU HAVE VARIABLES NAMES mean_est, mest_samples, Pnew

% compute mean, std, and norm for each sample
mean_samples  = zeros(Pnew,1);
std_samples   = zeros(Pnew,1);
norm_samples  = zeros(Pnew,1);
for ii=1:P
    m = mest_samples(:,ii);
    mdev = m - mean_est;    % note: remove mean of ALL samples
    %mean_samples(ii)  = 
    %std_samples(ii)   = 
    %norm_samples(ii)  = 
end

%==========================================================================
