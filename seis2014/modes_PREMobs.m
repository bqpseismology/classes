%
% modes_PREMobs.m
%
% Preparatory script for sumatra_modes_template.m
%

clear, close all, clc

ddir = './data/';

% spheroidal modes
dfile = [ddir 'prem_Smodes.txt'];
[nobs,~,lobs,T,Tstd] = textread(dfile,'%f%s%f%f%f','headerlines',6);
fmhz = 1./T*1e3;
disp('spheroidal mode observations (measured from seismograms):');
for ii=1:length(nobs)
   disp(sprintf('n = %i, l = %2i, T = %8.2f +/- %.2f s, f = %6.3f mHz',...
       nobs(ii),lobs(ii),T(ii),Tstd(ii),fmhz(ii))); 
end 

% toroidal modes
dfile = [ddir 'prem_Tmodes.txt'];
[nobsx,~,lobsx,Tx,Tstdx] = textread(dfile,'%f%s%f%f%f','headerlines',6);
for ii=1:length(nobs)
   disp(sprintf('n = %i, l = %2i, T = %8.2f +/- %.2f s',nobs(ii),lobs(ii),T(ii),Tstd(ii))); 
end
fmhzx = 1./Tx*1e3;

% plot
% note: 0.2-1.0 mHz is the range of Park et al. (2005), Figure 1
% note: 2S1 had not been observed in 1980
xmin = -0.5; xmax = 10.5; df=0.1; fsize=14;

% toroidal modes (same as in the modes HW)
figure; hold on;
plot([xmin xmax],0.2*[1 1],'k--');
plot([xmin xmax],1.0*[1 1],'k--');
scatter(lobsx,fmhzx,12^2,nobsx,'filled');
colorbar;
xlabel('degree, l'); ylabel('frequency, mHz');
title({'toroidal modes for PREM, colored by n',...
    'open circles = spheroidal modes for PREM'});
axis([xmin xmax 0 4]); caxis([0 3]); grid on;
for ii=1:length(nobsx)
    text(lobsx(ii),fmhzx(ii)+df,sprintf('%iT%i',nobsx(ii),lobsx(ii)),'fontsize',fsize);
end
scatter(lobs,fmhz,12^2,nobs,'ko');
fontsize(14);

% spheroidal modes
figure; hold on;
plot([xmin xmax],0.2*[1 1],'k--');
plot([xmin xmax],1.0*[1 1],'k--');
scatter(lobs,fmhz,12^2,nobs,'filled');
colorbar;
xlabel('degree, l'); ylabel('frequency, mHz');
title({'spheroidal modes for PREM, colored by n',...
    'open circles = toroidal modes for PREM'});
axis([xmin xmax 0 4]); caxis([0 3]); grid on;
for ii=1:length(nobs)
    text(lobs(ii),fmhz(ii)+df,sprintf('%iS%i',nobs(ii),lobs(ii)),'fontsize',fsize);
end
scatter(lobsx,fmhzx,12^2,nobsx,'ko');
fontsize(14);
