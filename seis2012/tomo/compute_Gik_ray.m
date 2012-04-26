%
% compute_Gik_ray.m
%
% Template script for computing an element of the partial derivatives matrix, Gik.
%

format short
format compact
close all
clear

%-----------------------------

ax1 = [-120.157113 -114.809623 32 36.364429];
lonmin = ax1(1); lonmax = ax1(2);
latmin = ax1(3); latmax = ax1(4);

%-----------------------------
% LOAD DATA

% load sources
[slon,slat,sinds] = textread('events_lonlat.dat','%f%f%f','headerlines',1);
nsrc = length(slat);

% load receivers
[rlon,rlat,rinds] = textread('recs_lonlat.dat','%f%f%f','headerlines',1);
nrec = length(rlat);

% load spline centers
[clon,clat] = textread('con_lonlat_q08.dat','%f%f','headerlines',0);
nspline = length(clat);

figure; hold on;
plot(clon,clat,'.');
text(clon,clat,num2str([1:nspline]'),'fontsize',6); % only seen when printed/saved
axis equal, axis(ax1);
xlabel(' Longitude'); ylabel(' Latitude');
title(' Center-points of spherical spline basis functions');
%fontsize(10); orient tall, wysiwyg

% REFERENCE HOMOGENEOUS PHASE VELOCITY
c0 = 3500;      % m/s

%-----------------------------
% compute ray paths (great circles)

% spline evaluations
opts = {1};
q = 8;

% number of points along each ray path
nump = 1000;

% test the ordering scheme for the rows of G
if (1==1)
    %i = 0;
    disp('     i   isrc  irec ');
    for isrc = 1:nsrc
        for irec = 1:nrec
            %i = i+1;
            i = (isrc-1)*nrec + irec;
            disp(sprintf('%6i%6i%6i',i,isrc,irec));
        end
    end
end

%==========================================================================
% compute design matrix






% save design matrix
%save([outdir 'Amat_ray'],'Amat_ray');

%==========================================================================

