% xy_gridsearch.m

clear, close all, clc

fsize = 12;
msize = 18;
ddir = './';

%----------------------

read_data;
whos obs_rngchg
[xvec,yvec] = plot_model(obs_rngchg);

% grid for 2D search
x0 = [19:0.2:22];
y0 = [21:0.2:23];
z0 = 2.58;
v0 = 0.0034;
iokay = find(isnan(obs_rngchg)==0);

misfit = NaN(length(x0),length(y0));
for k = 1:length(x0)
    x0(k)
    for l = 1:length(y0)
        syn_rngchg = mogi2insar(x0(k),y0(l),z0,v0,0,0);  
        misfit(k,l) = sum((obs_rngchg(iokay) - syn_rngchg(iokay)).^2);
    end
end

[indx indy] = find(misfit == min(misfit(:)));

disp(['Source X coordinate: ' num2str(x0(indx))]);
disp(['Source Y coordinate: ' num2str(y0(indy))]);

% plot cross section of misfit function
figure;
imagesc(x0,y0,misfit');
set(gca,'ydir','normal'); hold on;
plot(x0(indx),y0(indy),'kp','markersize',msize,'markerfacecolor','w');
colorbar
set(gca,'FontSize',12);
xlabel('Easting [km]','FontSize',14);
ylabel('Northing [km]','FontSize',14);
title('Misfit Function','FontSize',14);
axis equal;     % WARNING: only for comparing x vs y
axis tight;

% plot solution - no mask
mogi2insar(x0(indx),y0(indy),z0,v0,1,0);  
subplot(2,1,1); hold on; plot(x0(indx),y0(indy),'kp','markersize',msize,'markerfacecolor','w');
subplot(2,1,2); hold on; plot(x0(indx),y0(indy),'kp','markersize',msize,'markerfacecolor','w');

% plot solution - with mask
syn_rngchg = mogi2insar(x0(indx),y0(indy),z0,v0,1,1);  
subplot(2,1,1); hold on; plot(x0(indx),y0(indy),'kp','markersize',msize,'markerfacecolor','w');
subplot(2,1,2); hold on; plot(x0(indx),y0(indy),'kp','markersize',msize,'markerfacecolor','w');

% plot residual field
residuals = obs_rngchg - syn_rngchg;
plot_model(residuals);

%==========================================
% implement V-zs search here



%==========================================
