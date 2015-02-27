% Generates n copies of a complex 2-D Gaussian random field 
% from its spectral distribution matrix

% Input: 
% k = vector of k values
% m = vector of m values
% C = spectral distribution matrix. Assumed to have been generated
%     from meshgrid(K,M), so that it is indexed as C(m,k)
% n = number of copies to generate

% Output
% phi = complex Gaussian random field with covariance function 
%       mhifft(C)


function phi = grf2(k, m, C, n)

Nx = max(size(k));
Ny = max(size(m));
dk = k(2) - k(1);
dm = m(2) - m(1);
Periodx = 2*pi/dk;
Periody = 2*pi/dm;

Cmtx = repmat(C,[1,1,n]);  % Generate copies of the spectral matrix
A = randn([Ny,Nx,n]);   % N(0,1) random variables
B = randn([Ny,Nx,n]);

% Generate random field 
phi = sqrt(Periodx*Periody*Cmtx/2).*(A + 1i*B);


