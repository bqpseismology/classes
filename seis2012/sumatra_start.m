%
% sumatra_start.m
% Carl Tape, GEOS 694, Applied Seismology, March 2012
%
% Introduction to data analysis and the 2004 Mw 9.2 Sumatra-Andaman earthquake
%
% Adapted from run_getwaveform.m and getwaveform_input.m
%
%

clear
clc
close all

idatabase = 7;

% Sumatra parameters from GCMT catalog
[otime,tshift,hdur,lat,lon,dep,M,M0,Mw,eid] = readCMT_all;
imatch = find(and(Mw >= 8.9,year(otime)==2004));
datestr(otime(imatch))
otime_pde = datenum(2004,12,26,0,58,50);
otime_cmt = otime(imatch);
originTime = otime_pde;     % note: NOT centroid time
elat = lat(imatch);
elon = lon(imatch);
edep_km = dep(imatch);
eid = eid(imatch);
mag = Mw(imatch);

% note: BH_ also available
%chan = {'LHZ'};
chan = {'LHZ','LHE','LHN'};
%chan = {'BHZ'};
%chan = {'BHZ','BHE','BHN'};

duration_s = 3*3600;
tshift = 0.5*3600;
cutoff = [];
samplerate = [];
axb = [];

startTime = originTime - tshift/86400;
endTime   = originTime + duration_s/86400;
dur_dy = endTime-startTime;
fprintf('startTime is %s\n',datestr(startTime,31));
fprintf('total length of time requested: %.2f s (= %.2f min = %.2f hours)\n',...
    dur_dy*86400,dur_dy*3600,dur_dy*24);

% additional user parameters
sacdir = [];
iint = 0;            % integrate waveforms: =1 for displacement, =0 for velocity
iprocess = 1;

% get waveform object
[w,s,site,sitechan] = getwaveform(idatabase,startTime,endTime,chan,iint,iprocess,cutoff,samplerate,axb,sacdir,originTime,elat,elon,edep_km,mag,eid);
   
% plot record section
% default record section plotting parameters (only relevant if irs = 1)
% note: these can be over-ridden by specifications in getwaveform_input.m
isort = 2;
iabs = 1;
T1 = []; T2 = []; 
tshift = 0;
pmax = 50;
iintp = 0;
inorm = 1;
nfac = 1;
azcen = [];
iunit = 2;
tlims = [];

plotw_rs(w,isort,iabs,tshift,T1,T2,pmax,iintp,inorm,tlims,nfac,azcen,iunit);

%--------------------------------------------------------------------------

sta = get(w,'station');
imatch = find(strcmp('ESK',sta));
w0 = w(imatch);
w0 = w0([1 3 2]);   % order Z-N-E
get(w0,'channel')

figure; plot(w0(1));

figure; plot(integrate(detrend(w0(1))));

w1 = integrate(detrend(w0));

% threecomp object
TC = threecomp(w1);
figure; plot(TC);

TCrot = rotate(TC);
figure; plot(TCrot);

%==========================================================================
