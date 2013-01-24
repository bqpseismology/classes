%
% least_squares.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% This program introduces the least squares method for the example of
% fitting a line (i.e., a model with two parameters) to a set of scattered
% data.  We show three options for solving the problem, each which gives
% the same result. No regularization is needed (or used) in this example.
%

clear
close all
%format long
%format compact

%---------------------------

% USER PARAMETERS (CHANGE THESE)
n = 100;                % number of observations
sigvec = [0 0.3];       % standard deviations of added errors

% TARGET model vector (y-intercept, slope)
mtar = [2.1 -0.5]';
m = length(mtar);

% compute design matrix
x = linspace(-2,2,n)';  % input x-values
G = [ones(n,1) x];      % n by m design matrix
N = inv(G'*G)*G';       % m by n data resolution matrix
                        % (notice the matlab comment associated with inv)
                         
% display dimensions of these variables
whos

nsig = length(sigvec);
for kk = 1:nsig

    % generate errors -- these will be different for each run
    sigma = sigvec(kk);
    e = sigma * randn(n,1); % normally distributed random numbers
                            
    % generate target 'data'
    d = G*mtar + e;
    
    % optional: add one big anomaly
    %d(1) = d(1) + 1000*sigma;

    % SOLVE: compute least squares solution, estimates, and estimated variance.
    % (We show three options for mest, each with the same result.)
    mest = G\d              % estimated model
    %mest = N*d;
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

% search range, measured by the distance from the target model
m1_ran = 1;
m2_ran = 1;

% number controlling the number of gridpoints in model space
nx = 100;
%nx = 10;   % for gradient plot

% generate grid for the model space
m1_vec = linspace(mtar(1)-m1_ran, mtar(1)+m1_ran, nx);
m2_vec = linspace(mtar(2)-m2_ran, mtar(2)+m2_ran, nx);
[X,Y] = meshgrid(m1_vec,m2_vec);
[a,b] = size(X);
m1 = reshape(X,a*b,1);
m2 = reshape(Y,a*b,1);

% compute misfit function (and gradient)
G = [ones(n,1) x];              % design matrix
RSS = zeros(n,1);               % initialize misfit function

for kk=1:a*b
    mtry = [m1(kk) m2(kk)]';    % a sample from model space
    dtry = G*mtry;              % predictions from the model
    res = d - dtry;             % residuals between data and predictions
    RSS(kk) = sum(res.*res);    % residual sum of squares
    
    % COMPUTE GRADIENT HERE

end
Z = reshape(RSS,a,b);           % reshape for plotting

% plot the misfit function
nc = 30;    % number of contours to plot
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

% PLOT GRADIENT HERE WITH quiver COMMAND


%==========================================================================
