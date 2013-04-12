
clear, close all, clc

read_data;
x0 = [19:0.2:22];
y0 = [21:0.2:23];
z0 = 2.581;
v0 = 0.0033911;
iokay = find(isnan(obs_rngchg)==0);

misfit = NaN(length(x0),length(y0));
for k = 1:length(x0)
    x0(k)
    for l = 1:length(y0)
        modelout = mogi2insar(x0(k),y0(l),z0,v0,0);  
        misfit(k,l) = sum((obs_rngchg(iokay) - modelout(iokay)).^2);
    end
end

[indx indy] = find(misfit == min(misfit(:)));

disp(['Source X coordinate: ',num2str(x0(indx))]);
disp(['Source Y coordinate: ',num2str(y0(indy))]);

% plot cross section of misfit function
figure;
imagesc(x0,y0,misfit'); hold on;
colorbar
set(gca,'FontSize',12);
xlabel('Easting [km]','FontSize',14);
ylabel('Northing [km]','FontSize',14);
title('Misfit Function','FontSize',14);
plot(x0(indx),y0(indy),'kp','markersize',24,'markerfacecolor','w');
axis equal;     % WARNING: only for comparing x vs y
axis tight;

% plot solution
modelout = mogi2insar(x0(indx),y0(indy),z0,v0,1);  
subplot(2,1,1); hold on; plot(x0(indx),y0(indy),'kp','markersize',24,'markerfacecolor','w');
subplot(2,1,2); hold on; plot(x0(indx),y0(indy),'kp','markersize',24,'markerfacecolor','w');

% plot residual field
residuals = (obs_rngchg - (modelout));
plot_model(residuals);
