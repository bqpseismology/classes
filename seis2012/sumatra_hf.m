% sumatra_hf.m
% 
% Template script for analyzing the direct arrival waveforms from Sumatra:
%     channel BHZ, duration of up to 2 hours
%
% Adapted from run_getwaveform.m
%
%

clear
close all
clc

% CHANGE THIS TO YOUR BASE GEOTOOLS DIRECTORY
% add path (in principle, this only needs to be executed once)
gdir = '/home/carltape/GEOTOOLS/';
addpath([gdir 'matlab_util/']);

spdy = 86400;   % seconds per day

if 0==1
    % quick plot example
    startTime = 7.323055408564815e+05;
    endTime = 7.323175408564815e+05;
    ds = datasource('antelope','/home/admin/databases/SUMATRA/data/sumatra');
    scnl = scnlobject('CAN','LHZ','G','');
    w = waveform(ds,scnl,startTime,endTime);
    figure; plot(w,'xunit','h');
end

% extract the full database of BHZ waveforms
otime_pde = datenum(2004,12,26,0,58,50);
originTime = otime_pde;
startTime = originTime - 10;
endTime   = originTime + 10;

% the advantage of using getwaveform.m is that it will add all kinds of
% headers to the waveform objects, such as station azimuth, distance, etc
idatabase = 7;
channel = {'BHZ'};
iint = 0;
iprocess = 1;           % calibration applied (nm/s)
cutoff = []; samplerate = []; axb = []; sacdir = [];
% source parameters from GCMT catalog
elat = 3.09;
elon = 94.26;
edep_km = 28.6;
mag = 8.9974;
eid = 'M122604A';
w = getwaveform(idatabase,startTime,endTime,channel,iint,...
    iprocess,cutoff,samplerate,axb,sacdir,originTime,elat,elon,edep_km,mag,eid);
nw = length(w);

disp('here is a list of the waveforms you have:');
for ii=1:nw
   disp(sprintf('%3i %7s %3s %6s %10s',ii,get(w(ii),'channel'),get(w(ii),'KNETWK'),...
       get(w(ii),'station'),get(w(ii),'units')));
end

% save a copy to avoid rerunning
w0 = w;

% pick a subset of waveforms
%ipick = [1:nw];                    % default
ipick = [21 22 107 67 68 165];     % CHANGE THIS     
w = w0(ipick);

% PLOTTING PARAMETERS FOR plotw_rs.m (CHANGE THESE AS NEEDED)
isort = 2;      % =1 by azimuth, =2 by distance
iabs = 0;
tshift = [];
tmark = [];
pmax = 40;
iintp = 0;
inorm = 1;
tlims = [];
nfac = 1;
azcen = [];
iunit = 2;
imap = 0;

% plot record section
T1 = [];
T2 = [];
plotw_rs(w,isort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azcen,iunit,imap);

% plot map
[sta,rlat,rlon,elat,elon] = getm(w,'station','STLA','STLO','EVLA','EVLO');
plot_event_station(elat(1),elon(1),rlat,rlon,sta);

% SOME EXAMPLES OF USING THE PLOTTING COMMANDS

% example of cutting a record
w([4 6]) = [];
plotw_rs(w,isort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azcen,iunit,imap);

% example of applying a time shift
% note this is in the order of listed stations (NOT as ordered in the record section)
get(w,'station')
tshift = [1186 1250 845 1440];
plotw_rs(w,isort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azcen,iunit,imap);

% example of resetting plotting range
tlims = [-200 800];
plotw_rs(w,isort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azcen,iunit,imap);

% START YOUR ANALYSIS HERE


%==========================================================================