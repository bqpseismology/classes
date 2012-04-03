function f = fftvec(t)
% FFTVEC provides frequency vector for Matlab's fft convention
% 
% This is the bare-bones FFT example for teaching purposes in the sense
% that there is no padding with zeros, no odd-number time series, no
% all-positive frequencies, and no use of fftshift.
%
% EXAMPLE:
%   t = [2:0.4:8.0]'; n=length(t); f = fftvec(t); ix = 4; f([ix+1 n-ix+1])
%   dt = t(2)-t(1); df = 1/(2*dt)/(n/2) = 1/(dt*n)
%

n = length(t);
if mod(n,2)==1
   error('time series must have an even number of points');
end

dt = t(2) - t(1);       % sample rate
fNyq = 1/(2*dt);        % Nyquist frequency

% first half of frequency vector (note: first entry is f=0)
f1 = linspace(0, fNyq, n/2+1)';

% full frequency vector
f = [f1(1:n/2) ; -f1(end:-1:2)];

%==========================================================================
