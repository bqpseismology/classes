%
% tomography_hw.m
% Applied Seismology, April 2012, Carl Tape
% 
% Starting script for homework on seismic tomography.
%

clear
close all
format short
format compact

colors;

ax1 = [-121 -114 31 37];        % lon-lat plotting dimensions

%==========================================================================
% LOAD DATA

% load sources
[slon,slat,sind] = textread('events_lonlat.dat','%f%f%f','headerlines',1);
nsrc = length(slat);

% load receivers
[rlon,rlat,rind] = textread('recs_lonlat.dat','%f%f%f','headerlines',1);
nrec = length(rlat);

% load spline centers
q = 8;
[qlon,qlat] = textread('con_lonlat_q08.dat','%f%f','headerlines',0);
nspline = length(qlat);

%==========================================================================
% lon-lat gridpoints for plotting

numx = 100;
[lonplot,latplot] = gridvec(ax1(1),ax1(2),numx,ax1(3),ax1(4));
nplot = length(lonplot);

% Compute design matrix for expanding a function in terms of splines;
% this is needed to view the tomographic models that we generate at the end.
B = zeros(nplot,nspline);
for ii=1:nspline
    ff = spline_vals(qlon(ii),qlat(ii),q,lonplot,latplot,{1});
    B(:,ii) = ff(:);
end

%==========================================================================
% INVERSE PROBLEM



%==========================================================================
% PLOTTING THE SOLUTION(S)

% values from GMT 'seis' color palette (type 'colormap(seis)')
seis = [170	0	0;  206	0	0;  243	0	0;  255	24	0;  255	60	0;  255	97	0;
        255	133	0; 255	170	0;  255	206	0;  255	243	0;  255	255	0; 255	255	0;
        231	255	4;  161	255	17;  90	255	30;  51	249	64;  13	242	99;  0	194	152;
        0	125	214;  0	68	248;  0	34	226];
seis = cpt2cmap(seis);



%==========================================================================
