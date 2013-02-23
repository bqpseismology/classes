%
% forward_epicenter.m
% Carl Tape
% 18-March-2010
%
% This forward problem describes the straight-line ray path travel time
% from a source to a receiver in a homogeneous medium with unknown
% velocity.
%
% NOTE: If there are constants in the Matlab in-line functions, then these
% constants are assigned at the time the function is INITIALIZED, not at
% the time the function is called.
%
% CALLED BY : optimization.m
% CALLS     : plot_epicenters
%

%=========================================
% FORWARD PROBLEM    

% labels for parameters
mlabs = {'ts','xs','ys','v = ln(V/V0)'};
ulabs = {'s','km','km','none'};

% travel time computation (homogeneous velocity; straight ray paths)
% v is the logarithmic velocity, v = ln(V/V0), which is Cartesian
V0 = 1;
d2 = @(x,y,xr,yr) ( (xr-x)^2 + (yr-y)^2 );
d1 = @(x,y,xr,yr) ( sqrt(d2(x,y,xr,yr)) );
tt = @(t,x,y,v,xr,yr) (t + d1(x,y,xr,yr)/(V0*exp(v)) );

% receiver locations for uniform grid
xmin = 10; xmax = 80;
ymin = 20; ymax = 90;
xvec = linspace(xmin,xmax,4);
yvec = linspace(ymin,ymax,3);
[X,Y] = meshgrid(xvec,yvec);
xrec = X(:);
yrec = Y(:);

% computation of predictions
d = @(t,x,y,v) ([   tt(t,x,y,v,xrec(1),yrec(1))
                    tt(t,x,y,v,xrec(2),yrec(2))
                    tt(t,x,y,v,xrec(3),yrec(3))
                    tt(t,x,y,v,xrec(4),yrec(4))
                    tt(t,x,y,v,xrec(5),yrec(5))
                    tt(t,x,y,v,xrec(6),yrec(6))
                    tt(t,x,y,v,xrec(7),yrec(7))
                    tt(t,x,y,v,xrec(8),yrec(8))
                    tt(t,x,y,v,xrec(9),yrec(9))
                    tt(t,x,y,v,xrec(10),yrec(10))
                    tt(t,x,y,v,xrec(11),yrec(11))
                    tt(t,x,y,v,xrec(12),yrec(12))
                   ]);

% N x M matrix of partial derivatives (differentiate tt with respect to each parameter)
% evaluated at model m = (t,x,y,v)
G = @(t,x,y,v) ([
    1   -(d2(x,y,xrec(1),yrec(1)))^(-1/2)*(xrec(1)-x)/(V0*exp(v))   -(d2(x,y,xrec(1),yrec(1)))^(-1/2)*(yrec(1)-y)/(V0*exp(v))   -d1(x,y,xrec(1),yrec(1))/(V0*exp(v))
    1   -(d2(x,y,xrec(2),yrec(2)))^(-1/2)*(xrec(2)-x)/(V0*exp(v))   -(d2(x,y,xrec(2),yrec(2)))^(-1/2)*(yrec(2)-y)/(V0*exp(v))   -d1(x,y,xrec(2),yrec(2))/(V0*exp(v))
    1   -(d2(x,y,xrec(3),yrec(3)))^(-1/2)*(xrec(2)-x)/(V0*exp(v))   -(d2(x,y,xrec(3),yrec(3)))^(-1/2)*(yrec(2)-y)/(V0*exp(v))   -d1(x,y,xrec(3),yrec(3))/(V0*exp(v))
    1   -(d2(x,y,xrec(4),yrec(4)))^(-1/2)*(xrec(4)-x)/(V0*exp(v))   -(d2(x,y,xrec(4),yrec(4)))^(-1/2)*(yrec(4)-y)/(V0*exp(v))   -d1(x,y,xrec(4),yrec(4))/(V0*exp(v))
    1   -(d2(x,y,xrec(5),yrec(5)))^(-1/2)*(xrec(5)-x)/(V0*exp(v))   -(d2(x,y,xrec(5),yrec(5)))^(-1/2)*(yrec(5)-y)/(V0*exp(v))   -d1(x,y,xrec(5),yrec(5))/(V0*exp(v))
    1   -(d2(x,y,xrec(6),yrec(6)))^(-1/2)*(xrec(6)-x)/(V0*exp(v))   -(d2(x,y,xrec(6),yrec(6)))^(-1/2)*(yrec(6)-y)/(V0*exp(v))   -d1(x,y,xrec(6),yrec(6))/(V0*exp(v))
    1   -(d2(x,y,xrec(7),yrec(7)))^(-1/2)*(xrec(7)-x)/(V0*exp(v))   -(d2(x,y,xrec(7),yrec(7)))^(-1/2)*(yrec(7)-y)/(V0*exp(v))   -d1(x,y,xrec(7),yrec(7))/(V0*exp(v))
    1   -(d2(x,y,xrec(8),yrec(8)))^(-1/2)*(xrec(8)-x)/(V0*exp(v))   -(d2(x,y,xrec(8),yrec(8)))^(-1/2)*(yrec(8)-y)/(V0*exp(v))   -d1(x,y,xrec(8),yrec(8))/(V0*exp(v))
    1   -(d2(x,y,xrec(9),yrec(9)))^(-1/2)*(xrec(9)-x)/(V0*exp(v))   -(d2(x,y,xrec(9),yrec(9)))^(-1/2)*(yrec(9)-y)/(V0*exp(v))   -d1(x,y,xrec(9),yrec(9))/(V0*exp(v))
    1   -(d2(x,y,xrec(10),yrec(10)))^(-1/2)*(xrec(10)-x)/(V0*exp(v))   -(d2(x,y,xrec(10),yrec(10)))^(-1/2)*(yrec(10)-y)/(V0*exp(v))   -d1(x,y,xrec(10),yrec(10))/(V0*exp(v))
    1   -(d2(x,y,xrec(11),yrec(11)))^(-1/2)*(xrec(11)-x)/(V0*exp(v))   -(d2(x,y,xrec(11),yrec(11)))^(-1/2)*(yrec(11)-y)/(V0*exp(v))   -d1(x,y,xrec(11),yrec(11))/(V0*exp(v))
    1   -(d2(x,y,xrec(12),yrec(12)))^(-1/2)*(xrec(12)-x)/(V0*exp(v))   -(d2(x,y,xrec(12),yrec(12)))^(-1/2)*(yrec(12)-y)/(V0*exp(v))   -d1(x,y,xrec(12),yrec(12))/(V0*exp(v))
    ]);

% M x M matrix of second partial derivatives (only used in full Newton method)
% note: this contains the measurement index i
G2 = @(t,x,y,v,i)  ([
   0 0 0 0
   0  (d1(x,y,xrec(i),yrec(i)))^-3*(yrec(i)-y)^2/(V0*exp(v))            -(d1(x,y,xrec(i),yrec(i)))^-3*(xrec(i)-x)*(yrec(i)-y)/(V0*exp(v))   (d1(x,y,xrec(i),yrec(i)))^-1*(xrec(i)-x)/(V0*exp(v))
   0 -(d1(x,y,xrec(i),yrec(i)))^-3*(xrec(i)-x)*(yrec(i)-y)/(V0*exp(v))   (d1(x,y,xrec(i),yrec(i)))^-3*(xrec(i)-x)^2/(V0*exp(v))             (d1(x,y,xrec(i),yrec(i)))^-1*(yrec(i)-y)/(V0*exp(v))
   0  (d1(x,y,xrec(i),yrec(i)))^-1*(xrec(i)-x)/(V0*exp(v))               (d1(x,y,xrec(i),yrec(i)))^-1*(yrec(i)-y)/(V0*exp(v))                d1(x,y,xrec(i),yrec(i))/(V0*exp(v))
]);

%---------------------------------------------
% PRIOR MODEL (MEAN MODEL) : ts, xs, ys, v

% prior model
V = 5;
v = log(V/V0);              % logarithmic velocity (Cartesian)
mprior = [ 16 35 45 v ]';

% prior model covariance matrix (diagonal)
sigma_prior = [0.5 10 10 0.2]';         % standard deviations
cprior0     = diag( sigma_prior.^2 );   % diagonal covariance matrix
if inormalization==1
    Cmfac = nparm;
else
    Cmfac = 1;
end
cprior      = Cmfac * cprior0;          % WITH NORMALIZATION FACTOR
icprior     = inv(cprior);              % WITH NORMALIZATION FACTOR
icprior0    = inv(cprior0);
Lprior      = chol(cprior0,'lower')';   % square-root (lower triangular)

% sample the prior model distribution using the square-root UNNORMALIZED covariance matrix
for ii=1:nsamples, randn_vecs_m(:,ii) = randn(nparm,1); end
cov_samples_m  = Lprior * randn_vecs_m;
mprior_samples = repmat(mprior,1,nsamples) + cov_samples_m;

% compute the norm of each model sample using the inverse NORMALIZED covariance matrix
norm2_mprior = zeros(nsamples,1);
for ii=1:nsamples
    dm = mprior_samples(:,ii) - mprior;
    norm2_mprior(ii) = dm' * icprior * dm;
end
%figure; plot(norm2_mprior,'.')

%---------------------------------------------
% INITIAL MODEL

% minitial is DIFFERENT FOR EACH RUN, or you can fix it for testing purposes
if irandom_initial_model == 1
    minitial = mprior_samples(:,1);      % first sample
else
    minitial = [
        15.3890
        46.5236
        40.1182
        1.7748
    ];
end

%---------------------------------------------
% TARGET MODEL

% mtarget is DIFFERENT FOR EACH RUN, or you can fix it for testing purposes
if irandom_target_model == 1  
    mtarget = mprior_samples(:,end);      % last sample
else
    mtarget = [
           16.1314
           21.2922
           46.2974
            2.0903
                ];
end

%---------------------------------------------

% TARGET DATA
dtarget = d(mtarget(1),mtarget(2),mtarget(3),mtarget(4));

% data covariance matrix
sigma_obs = 0.5 * ones(ndata,1);    % standard deviations
cobs0     = diag( sigma_obs.^2 );   % diagonal covariance matrix
if inormalization==1
    Cdfac = ndata;
else
    Cdfac = 1;
end
cobs      = Cdfac * cobs0;          % WITH NORMALIZATION FACTOR
icobs     = inv(cobs);              % WITH NORMALIZATION FACTOR
icobs0    = inv(cobs0);
Lcobs     = chol(cobs0,'lower')';   % square-root (lower triangular)

% sample the data distribution using the square-root UNNORMALIZED covariance matrix
for ii=1:nsamples, randn_vecs_d(:,ii) = randn(ndata,1); end
cov_samples_d = Lcobs * randn_vecs_d;
dobs_samples  = repmat(dtarget,1,nsamples) + cov_samples_d;

% compute the norm of each data sample using the inverse NORMALIZED covariance matrix
norm2_dobs = zeros(nsamples,1);
for ii=1:nsamples
    dd = dobs_samples(:,ii) - dtarget;
    norm2_dobs(ii) = dd' * icobs * dd;
end
%figure; plot(norm2_dobs,'.')

% Pick the uncertainties for the target data by simply choosing a sample
%   from the realizations of the data covariance.
% NOTE: eobs is DIFFERENT FOR EACH RUN, or you can fix it for testing purposes
switch idata_errors
    case 0
        eobs = zeros(ndata,1);          % no errors
    case 1
        eobs = cov_samples_d(:,1);      % random: first sample
    case 2
        eobs = [                        % fixed
            -0.8689
            -0.4666
            -0.0516
            0.2411
            0.2825
            -0.2301
            0.1977
            0.3291
            1.0063
            0.5674
            0.1348
            0.4603
        ];
end

% "true" observations (includes added errors)
dobs = dtarget + eobs;

%---------------------------------------------
% PLOTS

% source-receiver geometry with ray paths
plot_epicenters([],mprior,minitial,mtarget,{xrec,yrec,1});
if iprint==1, print(gcf,'-depsc',[pdir 'srcrec_rays']); end

% source-receiver geometry with ray paths
plot_epicenters(mprior_samples,mprior,minitial,mtarget,{xrec,yrec,1});
if iprint==1, print(gcf,'-depsc',[pdir 'mprior_epi_rays']); end

% with prior samples (no ray paths)
opts = {xrec,yrec,0};
plot_epicenters(mprior_samples,mprior,minitial,mtarget,opts);
if iprint==1, print(gcf,'-depsc',[pdir 'mprior_epi']); end

%========================================================================
