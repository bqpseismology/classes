%
% CAN_response_template.m
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
spdy = 86400;   % seconds per day

%----------------------------------
% USER CHANGE THESE
iresponse = 1;  % Part 1 (lab_sumatra.pdf)
iwaveform = 0;  % Part 2 (Sumatra Homework, Part I)

% print figures to files
iprint = 0;
pdir = './';

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

%==========================================================================
% PART 1: EXAMINE THE INSTRUMENT RESPONSE
if iresponse==1

% specify frequencies
fmin = 1e-4;
fmax = 1e2;
f = logspace(log10(fmin),log10(fmax),100)';
omega = 2*pi*f;

aran = 10.^[-20 1.5];

% default option: get response from antelope database
% FIGURE 1
res0 = response_get_from_db(station,channel,startTime,f,dbname);
response_plot(res0); xlim([fmin fmax]); ylim(aran);
title('response_get_from_db.m','interpreter','none');
if iprint==1, print(gcf,'-depsc',sprintf('%sCAN_response_fig1',pdir)); end

% Here we create a response object by directly specifying the response file.
% For one file we manually removed the FIR filters in order to show that it
% would then match the response obtained from the converted-to-velocity
% sac pole-zero file.
% FIGURES 2 and 3
for kk=1:2
    if kk==1
        rfile0 = 'STRECKEISEN_STS1.5';
    else
        rfile0 = 'STRECKEISEN_STS1.5_noFIR';
    end
    rfile = [tdir rfile0];
    respObject = dbresponse(rfile);
    response.values = eval_response(respObject,omega);
    response.frequencies = f;
    response_plot(response);
    xlim([fmin fmax]); ylim(aran);
    title(rfile0,'interpreter','none');
    if iprint==1, print(gcf,'-depsc',sprintf('%sCAN_response_fig%i',pdir,kk+1)); end
end

% compare CAN response in antelope database (no FIR) with PZs from sac file
% FIGURE 4 (should match Figure 3)
ideriv = 1;     % velocity
[p,z,c,A0,k] = read_pzfile(pzfile,ideriv,1);   % velocity response
polezero.poles = p;
polezero.zeros = z;
polezero.normalization = A0;    % needed to match normalization
res = response_get_from_polezero(f,polezero);
response_plot(res); xlim([fmin fmax]); ylim(aran);
title(['sac pole-zero file: ' dlabs{ideriv+1}]);
if iprint==1, print(gcf,'-depsc',sprintf('%sCAN_response_fig4',pdir)); end

%-----------------------
% compare displacement, velocity, and acceleration spectra
% Note that all response plots are normalized such that the
% VELOCITY RESPONSE = 1 at some calibration period.
% FIGURE 5

% first: use read_pzfile.m (add one pole=0 per differentiation)
xf = 5;  % figure index
figure(xf); nr=3; nc=2;
for kk=1:3
    ideriv = kk-1;
    [p,z,c,A0,k] = read_pzfile(pzfile,ideriv);
    polezero.poles = p;
    polezero.zeros = z;
    polezero.normalization = c;    % note c, not A0
    res = response_get_from_polezero(f,polezero);
    % complex instrument response to either displacement, velocity, or acceleration
    Ix = res.values;
    % phase response
    subplot(nr,nc,2*kk-1); semilogx(f,angle(Ix)*deg); axis([fmin fmax -180 180]);
    xlabel('frequency, Hz'); ylabel('phase, deg');
    % amplitude response
    subplot(nr,nc,2*kk); loglog(f,abs(Ix)); xlim([fmin fmax]);
    xlabel('frequency, Hz'); ylabel('amplitude');
    title(sprintf('sac pole-zero file (%s)',dlabs{kk}));
end

% THIS IS ONLY A CHECK
% second: check this operation by dividing by i*omega for each operation
%         (here we plot in dashed red lines on Figure 10 to show an exact match)
% note: FFT[ f'(t) ] = i*w*FFT[ f(t) ]
[p,z,c,A0,k] = read_pzfile(pzfile,0);   % displacement
polezero.poles = p;
polezero.zeros = z;
polezero.normalization = c;
res = response_get_from_polezero(f,polezero);
Id = res.values;        % displacement
figure(xf); subplot(nr,nc,1); hold on; plot(f,angle(Id)*deg,'r--'); grid on;
figure(xf); subplot(nr,nc,2); hold on; plot(f,abs(Id),'r--'); axis(10.^[-4 2 -5 12]); grid on;
Iv = Id./(1i*omega);    % velocity
figure(xf); subplot(nr,nc,3); hold on; plot(f,angle(Iv)*deg,'r--'); grid on;
figure(xf); subplot(nr,nc,4); hold on; plot(f,abs(Iv),'r--'); axis(10.^[-4 2 -5 12]); grid on;
Ia = Id./(-omega.^2);   % acceleration
figure(xf); subplot(nr,nc,5); hold on; plot(f,angle(Ia)*deg,'r--'); grid on;
figure(xf); subplot(nr,nc,6); hold on; plot(f,abs(Ia),'r--'); axis(10.^[-4 2 -5 12]); grid on;
if iprint==1, orient tall; print(gcf,'-depsc',sprintf('%sCAN_response_fig%i',pdir,xf)); end

end  % iresponse

%==========================================================================
% PART 2: COMPUTE THE AMPLITUDE SPECTRUM OF ACCELERATION
if iwaveform==1
    
% load waveform
ds = datasource('antelope',dbname); 
scnl = scnlobject(station,channel,netwk,'');
w = waveform(ds,scnl,startTime,endTime);
wraw = remove_calib(w);
figure; plot(wraw); axis tight
if iprint==1, fontsize(14); print(gcf,'-depsc',sprintf('%sCAN_response_seis',pdir)); end

% get some info about the seismogram
[tstart,tdur,sps] = getm(w,'start','duration','freq');
dt = 1/sps;
disp(sprintf('duration is %.3f days = %.2f hours = %.2f min = %.3e s',...
    tdur,tdur*24,tdur*24*60,tdur*24*60*60));

% example of getting an absolute time from the seismogram
tpick = 3*1e5;      % based on the plot
datestr(tstart + tpick/spdy,31)

% for FFT (specifically for modes): demean and taper
w = demean(w); 
wd = get(w,'data');
nd = length(wd);
taper = tukeywin(nd,1);     % matlab taper function
wd = wd.*taper;
w = set(w,'DATA',wd);
sps = get(w,'freq');        % samples per second
dt = 1/sps;                 % time step

% compute FFT -- this can take 5-10 minutes
% note 1: this will attach the frequency version to the waveform object
% note 2: the file fftmat.m will be saved in your local directory,
%         so be sure to run Matlab from the same directory
fname = 'fftcan';
if ~exist([fname '.mat'],'file')
    tic, w = wf_fft.compute(w); toc
    save(fname,'w');
else
    load(fname);
end
f    = get(w,'fft_freq');
wAmp = get(w,'fft_amp');    % amplitude of H(w)
wPhs = get(w,'fft_phase');  % phase of H(w)

%--------------------------------------------------------------------------
% HOMEWORK/LAB EXERCISE STARTS HERE



%--------------------------------------------------------------------------

end  % iwaveform

%==========================================================================

