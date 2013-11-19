%
% lab_dispersion.m
% Carl Tape, Applied Seismology (GEOS 626)
% 
% Example of computing group dispersion and phase dispersion using two
% seismograms (real data).
%
% calls bandpass.m, markt.m, fftvec.m
%

clear
close all
clc
format compact

deg = 180/pi;
fsize = 8;

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
plot(ti,ynee,'r');
%plot(ti,yneeen,'b--',ti,-yneeen,'b--');    % envelope
xlim(xran);
xlabel('Time (s)'); ylabel('Amplitude'); title('Needles, LHZ');
%print(gcf,'-depsc',[pdir 'PAS_NEE_seis']);

%--------------------------------------------------------------------------

% CONSTRUCT FREQUENCY VECTOR
npt = nt;
f = fftvec(ti);
f = abs(f);      % we do not have negative frequencies here

Hp = fft(ypas);  %
Ap = abs(Hp);    % =sqrt(H.*conj(H)), where P=H.*conj(H) is the power spectral density
Hn = fft(ynee);
An = abs(Hn);

% explore these to see various details with the matlab FFT
if 0==1
    y = ypas; H = Hp; A = Ap;
    imax = npt/2+1;
    ip = 1:imax;
    
    % note: the first entry contains the INTEGRATED content of the signal
    disp('check the first entry of the FFT:');
    sum(ypas), mean(ypas)*npt, H(1)

    % check the difference between abs and power
    % if z = a + bi, then abs(z) = sqrt(z z*)
    norm( A.*A - H.*conj(H) ) / norm( A )

    % compare IFFT[FFT[y(t)]] with y(t)
    figure; plot(ti,y,'b',ti,ifft(H),'r--');

    % check the ordering of the complex entries of H
    figure; plot(real(H(2:imax-1)) - real(H(npt:-1:imax+1)),'.');
end

%-------------
% plot the spectrum for PAS and NEE
figure; hold on; 
plot(f,Ap,'b'); plot(f,An,'r');
legend('PAS','NEE');
xlabel('frequency (Hz)'); ylabel('spectral amplitude');
%print(gcf,'-depsc',[pdir 'PAS_NEE_spec']);

%==========================================================================
% CODE HERE FOR GROUP SPEED (use bandpass.m)


%==========================================================================
% CODE HERE FOR HARMONICS

Hp = fft(ypas);
Hn = fft(ynee);

whos ti ypas ynee f Hp Hn

tlims = [2600 2800];
numf = length(fvec);
for ii=1:numf
    ftar = fvec(ii);    % target frequency for harmonic
    
    % initialize fourier transforms
    Hp2 = complex(zeros(npt,1),0);
    Hn2 = complex(zeros(npt,1),0);
    
    % CODE HERE FOR HARMONICS
    

end

%==========================================================================
% CODE HERE FOR PHASE SPEED



%==========================================================================
