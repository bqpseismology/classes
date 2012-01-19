%
% hw1.m
% 
%
%

clc
clear
close all

% load GCMT catalog
idir = './data/';       % you might need to change the path to data
load([idir 'cmtall_sub']);
whos

% plot the catalog, colored by depth
[~,isort] = sort(dep,'descend');    % plot deepest events on top
msize = 4^2;                        % marker size
figure;
scatter(wrapTo360(lon(isort)),lat(isort),msize,dep(isort),'filled');
axis([0 360 -90 90]);
xlabel('Longitude'); ylabel('Latitutde');
caxis([0 600]); colorbar

figure;
plot(wrapTo360(lon(dep > 600)),lat(dep > 600),'.');
axis([0 360 -90 90]);
xlabel('Longitude'); ylabel('Latitutde');

%---------------------------

dmag = 0.1;
% note: seis2GR.m is not a built-in matlab command
[Ncum,N,Medges] = seis2GR(Mw,dmag);

