%
% function 
% CARL TAPE, 01-Dec-2009
% printed xxx
%
%
% calls xxx
% called by xxx
%

function [G1,G2] = chol_ogor(R,nx,ny,nsample)

n = length(R);         % matrix R is n by n

if n ~= nx*ny, error('R is not nx*ny by nx*ny'); end

% extract base Toeplitz matrices to make Rn (n x ny)
Rn = R(1:ny,:)';

% Cholesky decomposition
% CHOL(A) uses only the diagonal and upper triangle of A.
%     The lower triangle is assumed to be the (complex conjugate)
%     transpose of the upper triangle.  If A is positive definite, then
%     R = CHOL(A) produces an upper triangular R so that R'*R = A.
%     If A is not positive definite, an error message is printed.
[A,p] = chol(R);

% initialize output matrices
G1 = zeros(n,nsample);
G2 = zeros(n,nsample);

ii = 1;
%for ii=1:nsample

% Gaussian random vectors
if 0==1
    % use the same ones for each test
    if or(nx~=4,ny~=3) error('must have nx=4 and ny=3'); end
    phi_mat = [
        0.8631   -1.0639    0.0105   -0.4812
        0.6367   -1.2238    1.5133    1.0643
        -1.8474    0.4549   -1.7485   -0.0179 ];
else
    phi_mat = zeros(ny,nx);
    for jj=1:nx
        phi_mat(:,jj) = randn(ny,1);
    end
end
phi = phi_mat(:);       % n x 1

% use the Cholesky-derived square-root
G1(:,ii) = A*phi;

% initialize the output vector
g = zeros(n,1);

% start the algorithm
% --> "inv is slow and inaccurate"
R0 = Rn(1:ny,:)
R1 = Rn(ny+1:2*ny,:)
Qk = R0                     % Q0
Ck = chol(Qk,'lower')       % C0 (Ck*Ck' = Qk)
g(1:ny) = Ck*phi_mat(:,1)   % first ny entries of g

Bk = (R1*inv(R0))'          % B1 (ny x ny)
                            % Matlab: inv is slow and inaccurate
%Bk = diag(diag(R1))         % why does Bk = diag(R1)?
Bvk = Bk;

% loop over nx-1 Toeplitz block matrices
for k = 1:nx-1
    k
    %whos RvK Bk Jk Rkp1 Bkp1 g
    cind1 = 1:k*ny              % cumulative index
    ind0 = 1+(k-1)*ny:k*ny
    ind1 = k*ny+1:(k+1)*ny
    ind2 = (k+1)*ny+1:(k+2)*ny

    Rvk = Rn(cind1,:)

    R0, Rvk'*Bvk
    Qk = R0 - Rvk'*Bvk
    %Qk = Qk - Bvk(ind0,:)'*Qk*Bvk(ind0,:)

    Ck = chol(Qk,'lower')      % Qk is always ny x ny
    Jk = rot90(eye(ny)); jd = cell(1,k); [jd{:}] = deal(Jk); Jk = rot90(blkdiag(jd{:}))  % Jk = Ip

    g(ind1) = Bvk'*Jk*g(cind1) + Ck*phi_mat(:,k+1)

    if k < nx-1
        Rkp1 = Rn(ind2,:)
        Bkp1 = inv(Qk)*(Rkp1 - Rvk'*Jk*Bvk);  % Matlab: inv is slow and inaccurate
        Bvk = [(Bvk' - Bkp1*Bvk'*Jk)' ; Bkp1]
        % why are the blocks always (almost) diagonal?
        %ismall = find(abs(Bvk) < 1e-6); Bvk(ismall) = 0
    end
end

G2(:,ii) = g;

%end

%=========================================================================
% EXAMPLE

if 0==1
    clc, clear, close all, format long
    
    nx = 4; ny = 3; L = 1; sigmaM = 1; nsample = 6;
    n = nx*ny;
    stit = sprintf('n = %i, L = %i, sigma = %i',n,L,sigmaM);

    % generate mesh of points
    xmin = 1; xmax = nx; xinc = 1;
    ymin = 1; ymax = ny; yinc = 1;
    [x,y,d2,pA,pB] = xy2distance(xmin,xmax,xinc,ymin,ymax,yinc);

    % generate mesh of points
    D2 = reshape(d2,n,n); clear d2
    X = reshape(x,ny,nx);
    Y = reshape(y,ny,nx);

    R = sigmaM^2 * exp(-D2 / (2*L^2) );   % Gaussian covariance
    %R = sigmaM^2 * exp(-sqrt(D2) / L );   % exponential covariance
    G = chol_ogor(R,nx,ny,nsample);
    %disp(sprintf('n = %i, n^2 = %i, L = %i, sigmaM = %i, p = %i',n,n2,L,sigmaM,p ));

    figure; nr=3; nc=2;
    for ii=1:nsample
        subplot(nr,nc,ii);
        if 0==1
            pcolor(X, Y, reshape(G(:,ii),ny,nx) ); shading flat
        else
            imagesc(reshape(G(:,ii),ny,nx)); set(gca,'ydir','normal');
        end
        caxis(4*sigmaM*[-1 1]);
        %colorbar
        axis equal, axis([xmin-1 xmax+1 ymin-1 ymax+1]);
        %xlabel(' x '); ylabel(' y ');
        title(sprintf(' sample of Cm (%s)',stit));
    end
    orient tall, wysiwyg, fontsize(9)
end

% check that a single nx*ny Gaussian vector is the same as nx ny Gaussian
% vectors stacked together
if 0==1
    ny = 100; nx = 10;
    phi_mat = zeros(ny,nx);
    for jj=1:nx
        phi_mat(:,jj) = randn(ny,1);
    end
    phi = phi_mat(:);
    
    figure; nr=2; nc=2;
    subplot(nr,nc,1); plot_histo(phi,[-4:0.5:4]); ylim([0 0.3]); grid on;
    subplot(nr,nc,2); plot_histo(randn(nx*ny,1),[-4:0.5:4]); ylim([0 0.3]); grid on;
end

% generate J(n,p) of p. 17, Eq. 2
if 0==1
    n = 3;
    p = 4;
    Jp = rot90(eye(p));
    jd = cell(1,n);
    [jd{:}] = deal(Jp);
    Jd = rot90(blkdiag(jd{:}));
    
    for n=1:3
        Jp = rot90(eye(p)); jd = cell(1,n); [jd{:}] = deal(Jp); Jd = rot90(blkdiag(jd{:}))
    end
end

%=======================================================
