%
% Examples 2.1 and 2.2 
% from Parameter Estimation and Inverse Problems, 2nd edition, 2011
% by R. Aster, B. Borchers, C. Thurber
%

clear
close all
clc

% Load precomputed data
ddir = '/usr/local/matlab_toolboxes/aster/cd_5.2/Examples/chap2/ex_2_1_2/';
load([ddir 'data1.mat']);
t = data1(:,1);
y = data1(:,2);
sigma = data1(:,3);
N = length(t);
M = 3;

disp('displaying t,y,sigma')
[t , y , sigma ]

% build the parabolic system matrix
G = [ ones(N,1) t -0.5*t.^2 ];

% apply the weighting
yw = y ./ sigma;
Gw = G ./ [sigma sigma sigma];

% solve for the least-squares solution
disp('Least-squares solution')
m = inv(Gw'*Gw)*Gw'*yw

% get the covariance matrix
ginv = inv(Gw'*Gw)*Gw';
disp('Covariance matrix')
covm = ginv*ginv'

% get the 1.96-sigma (95%) conf intervals
disp('95% parameter confidence intervals (m-, mest, m+)')
dm = 1.96*sqrt(diag(covm));
[m-dm  m  m+dm]

% N-M degrees of freedom
dof = N-M;
disp(sprintf('Chi-square misfit for %i degrees of freedom (dof)',dof));
chi2 = norm((y - G*m)./sigma)^2

% Find the p-value for this data set
disp('chi-square p-value')
p = 1-chi2cdf(chi2,dof)

% Find the parameter correlations
s = sqrt(diag(covm))
disp('correlation matrix')
r = covm./(s*s');

% Plot the data and model predicted data
xx = min(t)-1 : 0.05 : max(t)+1;        % denser sampling of x values
mm = m(1) + m(2)*xx - 0.5*m(3)*xx.^2;

figure(1); hold on
plot(xx,mm,'k');
errorbar(t,y,sigma,'ko');
xlabel('Time (s)');
ylabel('Elevation (m)');
%print -deps2 c2fparabfig.eps

% Output covm and the eigenvalues/eigenvectors of covm.
disp('Covariance matrix for fitted parameters:')
covm
disp('Eigenvalues/eigenvectors of the covariance matrix:');
[u,lam] = eig(inv(covm))
disp('95% confidence ellipsoid semiaxis lengths:');
semi_axes = [sqrt(chi2inv(0.95,3)*(1./diag(lam)))]'

disp('95% confidence ellipsoid semiaxes:')
[semi_axes(1)*u(:,1), semi_axes(2)*u(:,2), semi_axes(3)*u(:,3)]

%-------------------------------

% Monte Carlo Section
y0 = G*m; 

nsamp = 1e3;

for nreal = 1:nsamp
  % Generate a trial data set of perturbed, weighted data
  ytrial = y0 + sigma.*randn(N,1);
  ywtrial = ytrial./sigma;
  mmc(nreal,:) = (Gw\ywtrial)';
  % Eq. 2.19
  chimc(nreal) =  norm( (G*mmc(nreal,:)' - ytrial) ./ sigma )^2;
end

% Plot the histogram of chi squared values
figure(2); nr=3; nc=1; xmax = 30; ymax = 0.14;

dbin = 0.5;
for ii=1:3
    subplot(nr,nc,ii);
    [B,Bplot] = plot_histo(chimc,[0:dbin:xmax],ii);
    sum(Bplot)*dbin  % check
    xlabel('\chi_{obs}^2');
    xlim([0 xmax]); if ii==3, ylim([0 ymax]); end
end

hold on;
xx = linspace(0,xmax,100);
dx = xx(2)-xx(1);
chitheo = chi2pdf(xx,dof);
sum(chitheo)*dx
plot(xx,chitheo,'r','linewidth',2);
title(sprintf('PDF for \\chi^2(\\nu=%i, x) PDF',dof))
%ylabel(sprintf('PDF for \\chi^2(\\nu=%i, x) PDF',dof))
%xlabel('x')
axis([0 xmax 0 ymax]);

% Plot the histograms of the model parameters
figure(3)
subplot(1,3,1); hist(mmc(:,1)); title('m_1 (m)')
subplot(1,3,2); hist(mmc(:,2)); title('m_2 (m/s)')
subplot(1,3,3); hist(mmc(:,3)); title('m_3 (m/s^2)')

% Plot the realizations of each pair of model parameters with the other
figure; nr=3; nc=2;
ax1 = [-50 50 85 110];
ax2 = [-50 50 7 12];
ax3 = [80 120 7 12];

subplot(nr,nc,1)
plot(mmc(:,1),mmc(:,2),'k.')
xlabel('m_1 (m)'); ylabel('m_2 (m/s)')
axis(ax1);

subplot(nr,nc,3)
plot(mmc(:,1),mmc(:,3),'k.')
xlabel('m_1 (m)'); ylabel('m_3 (m/s^2)')
axis(ax2);

subplot(nr,nc,5)
plot(mmc(:,2),mmc(:,3),'k.')
xlabel('m_2 (m/s)'); ylabel('m_3 (m/s^2)')
disp('Displaying Projections of 1000 Monte-Carlo models (fig 4)')
axis(ax3);

% Plot the 95% error ellipses for each pair of parameters

%generate a vector of angles from 0 to 2*pi
theta = (0:.01:2*pi)';
delta = sqrt(chi2inv(0.95,2));
%the radii in each direction from the center
r = zeros(length(theta),2);

%---------------

% compute the data for the m1, m2 ellipsoid.
C = covm((1:2),(1:2));
[u,lam] = eig(inv(C));
%calculate the x component of the ellipsoid for all angles
r(:,1) = (delta/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(delta/sqrt(lam(2,2)))*u(1,2)*sin(theta);
%calculate the y component of the ellipsoid for all angles
r(:,2) = (delta/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(delta/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m1, m2 ellipsoid
subplot(nr,nc,2)
plot(m(1)+r(:,1),m(2)+r(:,2),'k');
fill(m(1)+r(:,1),m(2)+r(:,2),'r');
axis(ax1);
xlabel('m_1 (m)'); ylabel('m_2 (m/s)');

% compute the data for the m1, m3 ellipsoid.
C = covm([1,3],[1,3]);
[u,lam] = eig(inv(C));
deltachisq = chi2inv(0.95,2);
delta = sqrt(deltachisq);
% calculate the x component of the ellipsoid for all angles
r(:,1) = (delta/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(delta/sqrt(lam(2,2)))*u(1,2)*sin(theta);
% calculate the y component of the ellipsoid for all angles
r(:,2) = (delta/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(delta/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m1, m3 ellipsoid
subplot(nr,nc,4)
plot(m(1)+r(:,1),m(3)+r(:,2),'k');
fill(m(1)+r(:,1),m(3)+r(:,2),'r');
axis(ax2);
xlabel('m_1 (m)'); ylabel('m_3 (m/s^2)');

% compute the data for the m2, m3 ellipsoid.
C = covm([2,3],[2,3]);
[u,lam] = eig(inv(C));
deltachisq = chi2inv(0.95,2);
delta = sqrt(deltachisq);
% calculate the x component of the ellipsoid for all angles
r(:,1) = (delta/sqrt(lam(1,1)))*u(1,1)*cos(theta)+(delta/sqrt(lam(2,2)))*u(1,2)*sin(theta);
% calculate the y component of the ellipsoid for all angles
r(:,2) = (delta/sqrt(lam(1,1)))*u(2,1)*cos(theta)+(delta/sqrt(lam(2,2)))*u(2,2)*sin(theta);

% plot the data for the m2, m3 ellipsoid
subplot(nr,nc,6)
plot(m(2)+r(:,1),m(3)+r(:,2),'k');
fill(m(2)+r(:,1),m(3)+r(:,2),'r');
axis(ax3);
xlabel('m_2 (m/s)'); ylabel('m_3 (m/s^2)');
%print -deps2 c2fellipseproj.eps

%==========================================================================

