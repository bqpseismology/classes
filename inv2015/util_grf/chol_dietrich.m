%
% function L = chol_dietrich(C,nx,ny)
% CARL TAPE, 01-Jan-2010
% printed xxx
%
% Fast Cholesky factorization of a block toeplitz matrix
%   C. R. Dietrich (1993a): 
%   C. R. Dietrich (1993b): 
%
% In fact, this algorithm appears to be somewhat slower than 
% Matlab's chol algorithm.
%
% KEY: THIS ALGORITHM ASSUMES THAT sigma = 1
%
% calls xxx
% called by xxx
%

function L = chol_dietrich(C,nx,ny)

% I define m = Nx and n = Ny, but I do not think this is what Dietrich meant;
% hence the re-assignment here.
m = ny;
n = nx;
k = m*n;
if length(C) ~= k, error('C is not m*n by m*n'); end
if m > n, error('must have m <= n'); end
disp(sprintf('(m,n,k) = (%i,%i,%i)',m,n,k));

disp(sprintf('C is %i by %i and comprises %i (%i by %i) Toeplitz blocks, each %i by %i',...
    k^2,k^2,nx*nx,nx,nx,ny,ny));

% set small values to zero, then estimate the bandwidth of the covariance matrix
% --> THIS DOES NOT SEEM TO BE HELPFUL
b = n;
if 0==1
    epscov2D = 1e-6;        % key parameter
    iwater = find(C < epscov2D);
    C(iwater) = 0;
    disp(sprintf('%i/%i entries of C are < %.2e',length(iwater),prod(size(C)),epscov2D));
    figure; spy(C);

    % estimate the bandwidth of the covariance matrix
    Ctop = C(:,1:m);
    b = n;
    for i=1:n
       Cblock = Ctop((i-1)*m+1:i*m,:);
       if ~any(Cblock(:))
           b = i;
           break
       end
    end
end
disp(sprintf('--> bandwidth is %i/%i',b,n));

% base matrix C1
%C1 = C(1:m,1:m); R1 = tril(C1); S1 = R1 - eye(m); C1,R1,S1

% initial decomposition
R = tril(C);
S = R - eye(k);
%norm( C - (R*R' - S*S') ) / norm(C)     % check

% initial U and V
U = R(:,1:m);
V = S(:,1:m);

tic

L = zeros(k,k);
for i = 0 : n-1
    nband = min([i+b n]);   % bandwidth
    %disp(sprintf('========== i = %i/%i ===========',i,n-1));
    %U,V
    for j1 = 1 : m
        for j2 = 1 : m
            rho = V(i*m+j1, j2) / U(i*m+j1, j1);
            alpha = 1/sqrt(1 - rho^2);
            %if ~isreal(alpha), disp(rho); end
            beta = rho*alpha;
            j = [i*m+1 : m*nband];
            %for j = i*m+1 : m*nband
                %disp(sprintf('(j,j1,j2) = (%i,%i,%i)',j,j1,j2));
                gamma = U(j,j1);
                U(j,j1) = alpha*gamma - beta*V(j,j2);
                V(j,j2) = -beta*gamma + alpha*V(j,j2);
            %end
        end
    end
    
    % assign U to columns of L
    L(:,i*m+1:(i+1)*m) = U;
    
    % the Zm operation
    U = [zeros(m,m) ; U(1:(n-1)*m,:) ];
    %U = [U(1:m,:) ; U(1:(n-1)*m,:) ];
end

toc

% ensure that L is lower triangular
L = tril(L);

% check if L is complex or not
% --> complex L is an indication of the instability of the algorithm
if ~isreal(L)
    icmplx = find( imag(L(:)) ~= 0);
    disp(sprintf('L has %i/%i (%.3f) complex entries',...
        length(icmplx),k^2,length(icmplx)/k^2));
    
    for i=1:k
        if ~isreal(L(:,i))
            icmplx = find( imag(L(:,i)) ~= 0);
            disp(sprintf('--> column %i/%i has %i complex entries',i,k,length(icmplx)));
        end
    end
    
    %L(icmplx) = 0;      % set complex entries to zero
    L = real(L);        % take the real part only
else
    disp('L is real');
end

%------------------------------------------------------

if 0==1
    clear, clc, close all
    nx = 4; ny = 3; Ls = 1;
    nx = 20; ny = 10; Ls = 2;
    nx = 20; ny = 10; Ls = 4;   % failure
    
    sigma = 0.05;
    k = nx*ny;
    stit = sprintf('k = %i, Ls = %i, sigma = %i',k,Ls,sigma);

    % generate mesh of points
    xmin = 1; xmax = nx; ymin = 1;
    [dx,ix0,iy0,iD2] = xy2distance(xmin,xmax,nx,ny); 
    
    % Gaussian covariance
    D2 = dx*iD2;
    C = sigma^2 * exp(-D2 / (2*Ls^2) );
    figure; imagesc(C); axis equal, axis tight
    
    L2 = sigma * chol_dietrich(C/sigma^2,nx,ny);
    figure; plot(L2(:)*sigma,'.')
    
    % compare with standard Cholesky algorithm
    tic, L1 = chol(C,'lower'); toc
    norm(L2), norm(L1), norm(L2 - L1) / norm(L1)
end

%=======================================================
