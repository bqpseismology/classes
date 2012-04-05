%
% response_CAN.m
% 
% This script shows some conventions for instrument response files
% associated with GISMO/Antelope and rdseed/sac.
% 
% calls numerous GISMO functions and also read_pzfile.m
%
% The example waveform is from CAN (Canberra, Australia) for the
% 2004 Mw 9.X Sumatra-Andaman earthquake.
%
% Carl Tape, 03/31/2012
%

clear
close all
clc

deg = 180/pi;

% set these to 1 or 0 for testing
iresponse = 1;
iwaveform = 0;

%----------------------------------

% test file is the 2004 Sumatra earthquake recorded at CAN.G, featured in
% Figure 1 of Park et al. 2005 (Science)
ddir = '/home/admin/databases/SUMATRA/data/';
tdir = [ddir 'CAN_test/'];
dbname = [ddir 'sumatra'];
station = 'CAN';
netwk = 'G';        % note: the network is not added to 'waveform' (get(w,'network'))
channel = 'LHZ';
stlab = [station '_' netwk '_' channel];
pzfile = [tdir 'SAC_PZs_G_CAN_LHZ__1989.153.00.00.00.0000_2006.344.02.60.60.99999'];

% waveform time interval
% note: start time needed to access response file
startTime0 = datenum(2004,12,25,12,58,50);
endTime0 = datenum(2005,1,4,12,58,50);
% pick times to get the entire waveform
startTime = startTime0 - 1;
endTime = endTime0 + 1;

dlabs = {'m to counts','m/s to counts','m/s^2 to counts'};

% PART 1; EXAMINE THE INSTRUMENT RESPONSE
if iresponse==1

% specify frequencies
fmin = 1e-4;
fmax = 1e2;
f = logspace(log10(fmin),log10(fmax),100)';
omega = 2*pi*f;

% default option: get response from antelope database
% FIGURE 1
res0 = response_get_from_db(station,channel,startTime,f,dbname);
response_plot(res0); xlim([fmin fmax]);
title('response_get_from_db.m','interpreter','none');

% Here we create a response object by directly specifying the response file.
% For one file we manually removed the FIR filters in order to show that it
% would then match the response obtained from the converted-to-velocity
% sac pole-zero file.
% FIGURES 2 and 3
for kk=1:2
    rfile0 = 'STRECKEISEN_STS1.5';
    if kk==1, rfile0 = 'STRECKEISEN_STS1.5_noFIR'; end
    rfile = [tdir rfile0];
    respObject = dbresponse(rfile);
    response.values = eval_response(respObject, omega);
    response.frequencies = f;
    response_plot(response);
    xlim([fmin fmax]); title(rfile0,'interpreter','none');
end

% compare CAN response in database with PZs from sac file
% FIGURE 4
ideriv = 1;
[p,z,c,A0,k] = read_pzfile(pzfile,ideriv,1);   % velocity response
polezero.poles = p;
polezero.zeros = z;
polezero.normalization = A0;    % needed to match normalization
res = response_get_from_polezero(f,polezero);
response_plot(res); xlim([fmin fmax]);
title(['sac pole-zero file: ' dlabs{ideriv+1}]);

%-----------------------
% compare displacement, velocity, and acceleration spectra
% response plots are normalized such that the VELOCITY RESPONSE = 1 at some
% calibration period
% FIGURE 10

% first use read_pzfile.m (add one pole=0 per differentiation)
figure(10); nr=3; nc=2;
for kk=1:3
    ideriv = kk-1;
    [p,z,c,A0,k] = read_pzfile(pzfile,ideriv);
    polezero.poles = p;
    polezero.zeros = z;
    polezero.normalization = c;    % note c, not A0
    res = response_get_from_polezero(f,polezero);
    H = res.values;
    subplot(nr,nc,2*kk-1); semilogx(f,angle(H)*deg); axis([fmin fmax -180 180]);
    xlabel('frequency, Hz'); ylabel('phase, deg');
    subplot(nr,nc,2*kk); loglog(f,abs(H)); xlim([fmin fmax]);
    xlabel('frequency, Hz'); ylabel('amplitude');
    title(sprintf('sac pole-zero file (%s)',dlabs{kk}));
end

% second: check this operation by dividing by i*omega for each operation
% note: FFT[ f'(t) ] = i*w*FFT[ f(t) ]
[p,z,c,A0,k] = read_pzfile(pzfile,0);   % displacement
polezero.poles = p;
polezero.zeros = z;
polezero.normalization = c;
res = response_get_from_polezero(f,polezero);
H = res.values;
P = H;
figure(10); subplot(nr,nc,1); hold on; plot(f,angle(P)*deg,'r--');
figure(10); subplot(nr,nc,2); hold on; plot(f,abs(P),'r--');
P = H./(1i*omega);
figure(10); subplot(nr,nc,3); hold on; plot(f,angle(P)*deg,'r--');
figure(10); subplot(nr,nc,4); hold on; plot(f,abs(P),'r--');
P = H./(-omega.^2);
figure(10); subplot(nr,nc,5); hold on; plot(f,angle(P)*deg,'r--');
figure(10); subplot(nr,nc,6); hold on; plot(f,abs(P),'r--');

end  % iresponse

%=====================================

% PART 2: COMPUTE THE AMPLITUDE SPECTRUM OF ACCELERATION
if iwaveform==1
    
% load waveform
% FIGURE 1
ds = datasource('antelope',dbname); 
scnl = scnlobject(station,channel ,netwk,'');
w = waveform(ds,scnl,startTime,endTime);
w = remove_calib(w);
figure; plot(w); axis tight
tstart = get(w,'start');

% example of getting an absolute time from the seismogram
tpick = 3*1e5;
datestr(tstart + tpick/86400,31)

% for FFT (specifically for modes): demean and taper
w = demean(w); 
wd = get(w,'data');
nd = length(wd);
taper = tukeywin(nd,1);     % matlab function
wd = wd.*taper;
w = set(w,'DATA',wd);
fNyq = get(w,'NYQ');

% compute FFT -- this can take awhile
% note 1: this will attach the frequency version to the waveform object
% note 2: the file fftmat.m will be saved in your local directory
%w = resample(w,'median',10);   % will speed up the FFT but not look as good
fname = 'fftcan';
if ~exist([fname '.mat'],'file')
    tic, w = wf_fft.compute(w); toc
    save(fname,'w');
else
    load(fname);
end
f    = get(w,'fft_freq');
wAmp = get(w,'fft_amp');
wPhs = get(w,'fft_phase');

%----------------------------------
% HOMEWORK/LAB EXERCISE STARTS HERE
% (1) describe what you 
% (1) plot the spectrum of the raw seismogram (use loglog)
% (2) construct the instrument response for acceleration
% (3) deconvolve the instrument response from the raw spectral seismogram
% (4) compare the spectra for the raw seismogram and instrument-deconvolved
%     seismogram over the frequency range [0.2 1.0] mHz



%----------------------------------

end  % iwaveform

%==========================================================================

