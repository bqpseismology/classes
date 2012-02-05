%
% least_squares.m
% Carl Tape, GEOS 694, Applied Seismology, Jan 2012
%
% This program introduces the least squares method for the example of
% fitting a line (i.e., a model with two parameters) to a set of scattered
% data.  We show three options for solving the problem, each which gives
% the same result. No regularization is needed (or used) in this example.
%

clear
close all
format long
format compact

%---------------------------

% TARGET model vector (y-intercept, slope)
mtar = [2.1 -0.5]';
m = length(mtar);

% compute design matrix
n = 100;                % USER CHANGE THIS: number of observations
x = linspace(-2,2,n)';  % input x-values
G = [ones(n,1) x];      % n by m design matrix
N = inv(G'*G)*G';       % m by n data resolution matrix
                        % (notice the matlab comment associated with inv)
                         
% display dimensions of these variables
whos

sigvec = [0 0.3];       % USER CHANGE THIS: error values

nsig = length(sigvec);
for kk = 1:nsig

    % generate errors -- these will be different for each run
    sigma = sigvec(kk);
    e = sigma * randn(n,1); % normally distributed random numbers
                            
    % generate target 'data'Jan
    d = G*mtar + e;
    
    % optional: add one big anomaly
    %d(1) = d(1) + 1000*sigma;

    % SOLVE: compute least squares solution, estimates, and estimated variance.
    % (We show three options for mest, each with the same result.)
    mest = N*d;             % estimated model
    %mest = G\d
    %mest = flipud(polyfit(x,d,1)')
    
    dest = G*mest;          % estimated predictions
    res = d - dest;         % residuals

    figure; msize = 10;
    stres = [' std(res) = ' sprintf('%.3f', std(res) )];

    subplot(2,1,1); hold on;
    plot(x,d,'.','markersize',msize);
    plot(x,dest,'r--','linewidth',2);
    xlabel(' x'); ylabel(' d');
    title({sprintf('Estimated model : m = (%.2f, %.2f)',mest(1),mest(2)), stres})
    grid on; axis equal;
    axis([min(x) max(x) min(G*mtar)-2*sigma max(G*mtar)+2*sigma]);

    subplot(2,2,3);
    plot(res,'.','markersize',msize); grid on; ylim([-1 1]);
    xlabel(' Observation index'); ylabel(' Residual, d - dest'); title(stres);

    subplot(2,2,4);
    edges = [-1.05:0.1:1.05]; [Nh,bin] = histc(res,edges);
    bar(edges,Nh,'histc'); xlim([min(edges) max(edges)]);
    xlabel(' Residual'); ylabel(' Number'); title([' Ntotal = ' num2str(n)]);

    %fontsize(11); orient tall, wysiwyg
end

break

%---------------------------
% generate a plot showing the RSS as a function of model space
% note: this is the 'brute-force' approach to solving this problem

% search range, measured by the distance from the target model
m1_ran = 1;
m2_ran = 1;

% generate a fine grid for the misfit function;
% use the target model for the center of the grid
npts = 100;
m1_vec = linspace(mtar(1)-m1_ran, mtar(1)+m1_ran, npts);
m2_vec = linspace(mtar(2)-m2_ran, mtar(2)+m2_ran, npts);
[X,Y] = meshgrid(m1_vec,m2_vec);
[a,b] = size(X);
m1 = reshape(X,a*b,1);
m2 = reshape(Y,a*b,1);

% compute misfit function
G = [ones(n,1) x];
RSS = zeros(n,1);
for kk=1:a*b
    mtry = [m1(kk) m2(kk)]';
    dtry = G*mtry;
    res = d - dtry;
    RSS(kk) = sum(res'*res);
end
Z = reshape(RSS,a,b);

% plot the misfit function
nc = 30;
figure; hold on;
contourf(X,Y,Z,nc); shading flat;
%scatter(m1,m2,6^2,RSS,'filled'); shading flat;
l1 = plot(mtar(1),mtar(2),'ws','markersize',10,'markerfacecolor','k');
l2 = plot(mest(1),mest(2),'wo','markersize',10,'markerfacecolor','r');
legend([l1,l2],'target model','estimated model');
axis equal, axis tight;
caxis([-1e-6 0.5*max(RSS)]); colorbar
xlabel(' m0, y-intercept');
ylabel(' m1, slope');
title(' Residual sum of squares');

% generate a coarse grid for the gradient
npts = 10;
m1_vec = linspace(mtar(1)-m1_ran, mtar(1)+m1_ran, npts);
m2_vec = linspace(mtar(2)-m2_ran, mtar(2)+m2_ran, npts);
[X,Y] = meshgrid(m1_vec,m2_vec);
[a,b] = size(X);
m1 = reshape(X,a*b,1);
m2 = reshape(Y,a*b,1);

% compute gradient gamma(m), then superimpose on F(m) plot


%==========================================================================
