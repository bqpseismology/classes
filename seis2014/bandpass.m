function y2 = bandpass(t,y,fa,fb)
% BANDPASS apply bandpass filter on a time series
%
% Code from Muller and MacDonald, Ice Ages and Astronomical Causes, p. 294
%

f1 = min([fa fb]);
f2 = max([fa fb]);
n = length(t);
dt = t(2) - t(1);
fNyq = 1/(2*dt);
y = y - mean(y);
ft = fft(y);
fre = linspace(0, 2*fNyq, n)';
[~, k1] = min(abs(fre - f1));
[~, k2] = min(abs(fre - f2));
k3 = n - k2 + 2;
k4 = n - k1 + 2;

ft(1:k1)  = zeros(k1, 1);
ft(k2:k3) = zeros(k3-k2+1, 1);
ft(k4:n)  = zeros(k1-1, 1);

y2 = real(ifft(ft));

%==========================================================================
