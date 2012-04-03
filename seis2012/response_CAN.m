%
% response_CAN.m
% 
% This script shows some conventions for instrument response files
% associated with GISMO/Antelope and rdseed/sac.
% 
% calls numerous GISMO functions and also read_pzfile.m
%
% The example waveform is from CAN (Canberra, Australia) for the 2004 Mw
% 9.X Sumatra-Andaman earthquake.
%
% Carl Tape, 03/31/2012
%

clear
close all
clc

deg = 180/pi;

iresponse = 0;
iwaveform = 1;

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

if iresponse==1

% specify frequencies
fmin = 1e-4;
fmax = 1e2;
f = logspace(log10(fmin),log10(fmax),100)';
omega = 2*pi*f;

% default option: get response from database
res0 = response_get_from_db(station,channel,startTime,f,dbname);
response_plot(res0); xlim(10.^[-4 0]);
title('response_get_from_db','interpreter','none');

% Here we create a response object by directly specifying the response file.
% For one file we manually removed the FIR filters in order to show that it
% would then match the response obtained from the converted-to-velocity
% sac pole-zero file.
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
[p,z,c,A0,k] = read_pzfile(pzfile,1);   % velocity response
polezero.poles = p;
polezero.zeros = z;
polezero.normalization = A0;    % needed to match normalization
res = response_get_from_polezero(f,polezero);
response_plot(res); xlim([fmin fmax]); title('sac PZ file');

%-----------------------
% compare displacement, velocity, and acceleration spectra
% response plots are normalized such that the VELOCITY RESPONSE = 1 at some
% calibration period

% first use read_pzfile.m (add one pole=0 per differentiation)
figure(10); nr=3; nc=2;
slabs = {'M to counts','M/S to counts','M/S^2 to counts'};
for kk=1:3
    ideriv = kk-1;
    [p,z,c,A0,k] = read_pzfile(pzfile,ideriv);
    polezero.poles = p;
    polezero.zeros = z;
    polezero.normalization = c;    % note c, not A0
    res = response_get_from_polezero(f,polezero);
    H = res.values;
    subplot(nr,nc,2*kk-1); semilogx(f,angle(H)*deg); xlim([f(1) f(end)]);
    xlabel('frequency, Hz'); ylabel('phase, deg');
    subplot(nr,nc,2*kk); loglog(f,abs(H)); xlim([f(1) f(end)]);
    xlabel('frequency, Hz'); ylabel('amplitude');
    title(sprintf('sac PZ file, ideriv=%i (%s)',ideriv,slabs{kk}));
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

% COMPUTE THE AMPLITUDE SPECTRUM OF ACCELERATION
if iwaveform==1
    
% load waveform
ds = datasource('antelope',dbname); 
scnl = scnlobject(station,channel ,netwk,'');
w = waveform(ds,scnl,startTime,endTime);
w = remove_calib(w);
figure; plot(w); 

% read sac file
sacval = rsac([tdir '2004.360.12.36.12.6559.G.CAN..LHZ.D.SAC']);
lh(sacval)
figure; plot(sacval(:,1),sacval(:,2));
title(['sac file for ' stlab],'interpreter','none');
xlabel('Times, s'); ylabel('Counts');

% check that they are identical
fwav = get(w,'data');
%ti = get(w,'timevector');   % in matlab days
ti = sacval(:,1);            % in seconds
fsac = sacval(:,2);
n = length(ti);
norm(fwav-fsac)             % check

% for FFT (specifically for modes): demean and taper
w = demean(w); 
wd = get(w,'data');
nd = length(wd);
taper = tukeywin(nd,1);
wd = wd.*taper;
w = set(w,'DATA',wd);
fNyq = get(w,'NYQ');

% compute FFT -- this can take awhile
% note: this will attach the frequency version to the waveform object
%w = resample(w,'median',10);
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

end  % iwaveform

%==========================================================================

