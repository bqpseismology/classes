%
% backus.m
%
% Implementation of the problem in Section 7.19 of Tarantola (2005),
%   "Using the Backus and Gilbert Method"
%

close all, clc, clear

zdep = [2.3 8.1 10.1];         % m
dobs = [2.0 8.1 10.2]';        % ms


G = @(z,i) ( and(z >= 0, z < zdep(i)) );

zmin = -1;
zmax = 15;
nz = 1000;
zplot = linspace(zmin,zmax,nz)';

figure; nr=3; nc=2;
ax0 = [zmin zmax -0.1 1.1];

subplot(nr,nc,1); hold on;
plot(zplot,G(zplot,1),'b');
plot(zplot,G(zplot,2),'r--');
plot(zplot,G(zplot,3),'k');
legend('G_1(z)','G_2(z)','G_3(z)');
axis(ax0);
xlabel('Depth, m'); ylabel('Slowness, s/km');
title('kernels of G');

S = diag(zdep);
S(1,2) = zdep(1);
S(1,3) = zdep(1);
S(2,3) = zdep(2);
S(2,1) = S(1,2);
S(3,1) = S(1,3);
S(3,2) = S(2,3);

Sinv = inv(S);

p = Sinv*dobs;
p = [-0.182 0.002 1.050]';

S, Sinv, p

nplot = [G(zplot,1) G(zplot,2) G(zplot,3)]*p;

subplot(nr,nc,2);
plot(zplot,nplot); axis(ax0);
xlabel('Depth, m'); ylabel('Slowness, s/km');

% resolving kernels
ax1 = [ax0(1:2) -0.05 0.6];
ztar = [1 5 10];
for k = 1:length(ztar)
    Rz = zeros(nz,1);
    for i=1:3
        for j=1:3
            Rz = Rz + G(ztar(k),i)*Sinv(i,j)*G(zplot,j);
        end
    end
    subplot(nr,nc,2+k);
    plot(zplot,Rz);
    axis(ax1);
    xlabel('z'', depth, m');
    ylabel(sprintf('R(z=%.1f m, z'')',ztar(k)));
end


