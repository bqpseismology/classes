function [f_r, rss, f_r_ss] = tsvd(g, X, r)
%TSVD truncated SVD regularization.
%
%    Given a vector g, a design matrix X, and a truncation parameter r, 
%           
%           [f_r, rss, f_r_ss] = tsvd(g, X, r) 
% 
%    returns the truncated SVD estimate of the vector f in the linear
%    regression model
% 
%            g = X*f + noise.
%
%    Also returned are the residual sum of squares rss and the sum of
%    squares f_r_ss of the elements of the truncated SVD estimate f_r
%    (the squared norm of f_r).
%
%    If r is a vector of truncation parameters, the i-th column
%    f_r(:,i) is the truncated SVD estimate for the truncation
%    parameter r(i); the i-th elements of rss and f_r_ss are the
%    associated residual sum of squares and estimate sum of squares.
%
%    Adapted from the TSVD routine in Per Christian Hansen's
%    Regularization Toolbox. 

% size of inputs
[n, p]      = size(X);
q           = min(n, p);
nr          = length(r);

% Possible choice of truncation parameter?
if ( min(r) < 0 || max(r) > q )
    error('Impossible truncation parameter r.')
end

% initialize outputs
f_r         = zeros(p, nr);
rss         = zeros(nr, 1); 
f_r_ss      = zeros(nr, 1);

% compute SVD of X
[U, S, V]   = svd(X, 0);  
s           = diag(S);      % vector of singular values

% "Fourier" coefficients (fc) in expansion of solution
% in terms of right singular vectors
beta        = U(:, 1:q)'*g;             % note data g
fc          = beta ./ s;      

% treat each truncation parameter separately
for j = 1:nr
	k         = r(j);                   % current truncation parameter
	f_r(:, j) = V(:, 1:k) * fc(1:k);    % truncated SVD estimated model vector
	f_r_ss(j) = sum(fc(1:k).^2);        % the squared norm of f_r
	rss(j)    = sum(beta(k+1:q).^2);    % norm(coefs*s_i NOT used)
end

% in overdetermined case, add rss of least-squares problem
if (n > p)
    rss = rss + sum((g - U(:, 1:q)*beta).^2);   % note data g
end

%==========================================================================
  