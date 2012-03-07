%
% disper.m
% Carl Tape, 06-March-2012
% Applied Seismology
% 
% Example of computing group dispersion and phase dispersion using two
% seismograms (real data).
%
% calls bandpass.m, markt.m
%

clear
close all
clc
format compact

deg = 180/pi;
fsize = 8;

iharm = 1;

%--------------------------------------------------------------------------

% target periods for measurements
Tvec = [20 30 40 50]';
fvec = 1./Tvec;
numf = length(fvec);

ax1 = [18 52 2.8 4.6];
delx = 331;             % distance between PAS and NEE (km)
    
% we are told that the phase velocity at each period must fall within these ranges
cran = [3.1 3.9; 3.0 4.3; 3.3 4.5; 3.3 4.5];

%--------------------------------------------------------------------------

% load data files
ww1 = 'pas.dat';
ww2 = 'nee.dat';
str1 = '/home/carltape/classes/caltech/kanamori_ge162/hw3/';
load([str1 ww1]);
load([str1 ww2]);

ti = nee(:,1);
dt = ti(2) - ti(1);
ynee = nee(:,2);
ypas = pas(:,2);
% for FFT an even number is extremely helpful, so chop off the first point
ti = ti(2:end);
ynee = ynee(2:end);
ypas = ypas(2:end);
nt = length(ti);

% use hilbert transform to compute envelope
ypasen = abs(hilbert(ypas));
yneeen = abs(hilbert(ynee));

figure; nr=2; nc=1;
xran = [ti(1) ti(end)];

subplot(nr,nc,1); hold on;
plot(ti,ypas,'b');
%plot(ti,ypasen,'r--',ti,-ypasen,'r--');    % envelope
xlabel('Time (s)'); ylabel('Amplitude'); title('Pasadena, LHZ');
xlim(xran);

subplot(nr,nc,2); hold on;
plot(ti,ynee,'b');
%plot(ti,yneeen,'r--',ti,-yneeen,'r--');    % envelope
xlim(xran);
xlabel('Time (s)'); ylabel('Amplitude'); title('Needles, LHZ');

%--------------------------------------------------------------------------

npt = nt;
fNyq = 1/(2*dt);                   % Nyquist frequency
f1 = linspace(0, fNyq, npt/2+1)';  % first half of frequency vector
                                   % note: first entry is f=0
f = [f1 ; f1(end-1:-1:2)];

H = fft(ypas);      % PAS
P = H.*conj(H);     % power spectral density
imax = npt/2+1;
ip = 1:imax;
figure; plot(f(ip),P(ip),'r');

% note: the first entry contains the INTEGRATED content of the signal
sum(ypas), mean(ypas)*npt, H(1)

% check the difference between abs and power
norm( abs(H).*abs(H) - H.*conj(H) ) / norm( H.*conj(H) )

% check
figure; hold on;
plot(ti,ypas,'b',ti,ifft(H),'r--');

% check the ordering
figure;
plot(real(H(2:imax-1)) - real(H(npt:-1:imax+1)),'.');

%-------------
% plot the spectrum for PAS
figure; nr=2; nc=2; fpmax = 0.1;
Z = abs(H);
ip = find(f <= fpmax);
ip = ip(2:ceil(length(ip)/2));  % one copy only; exclude f=0

subplot(nr,nc,1); plot(f,Z,'b');
xlabel('frequency (Hz)'); ylabel('spectral amplitude');
subplot(nr,nc,2); plot(f(ip),Z(ip),'b.-');
xlabel('frequency (Hz)'); ylabel('spectral amplitude');

subplot(nr,nc,3);
plot(log10(f(2:end)),log10(Z(2:end)),'b');
axis tight; ylim([4 7]);
xlabel('log10 frequency (Hz)'); ylabel('log10 amplitude');

subplot(nr,nc,4); 
plot(log10(f(ip)),log10(Z(ip)),'b.-');
axis tight; ylim([4 7]);
xlabel('log10 frequency (Hz)'); ylabel('log10 amplitude');
%-------------

% HARMONICS
Hp = fft(ypas);
Hn = fft(ynee);

whos ti ypas ynee f Hp Hn

tlims = [2600 2800];
numf = length(fvec);
for ii=1:numf
    ftar = fvec(ii);    % target frequency for harmonic

    % CODE HERE FOR HARMONICS
end

%==========================================================================
% CODE HERE FOR GROUP SPEED


%==========================================================================
% CODE HERE FOR HARMONICS

Hp = fft(ypas);
Hn = fft(ynee);

tlims = [2600 2800];
numf = length(fvec);
for ii=1:numf
    ftar = fvec(ii);    % target frequency for harmonic

    % CODE HERE FOR HARMONICS
    
end

%==========================================================================
% CODE HERE FOR PHASE SPEED



%==========================================================================
