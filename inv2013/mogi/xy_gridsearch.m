
clear, close all, clc

read_data;
x0 = [19:0.2:22];
y0 = [21:0.2:23];
z0 = 2.581;
v0 = 0.0033911;
temp = isnan(obs_rngchg);
a = find(temp == 0);

for k = 1:length(x0)
    x0(k)
    for l = 1:length(y0)
        modelout = mogi2insar(x0(k),  y0(l),   z0, v0,0);  
        misfit(k,l) = sum((obs_rngchg(a) - modelout(a)).^2);
    end
end

[indx indy] = find(misfit == min(misfit(:)));

disp(['Source X coordinate: ',num2str(x0(indx))]);
disp(['Source Y coordinate: ',num2str(y0(indy))]);

figure;imagesc(x0,y0,misfit');colorbar
set(gca,'FontSize',12); h = xlabel('Easting [km]');set(h,'FontSize',14);
h = ylabel('Northing [km]');set(h,'FontSize',14);
h = title('Misfit Function');set(h,'FontSize',14);
h = text(x0(indx),y0(indy),'*');set(h,'FontSize',24,'Color',[1 1 1]);
axis equal;     % WARNING: only for comparing x vs y
axis tight;

modelout = mogi2insar(x0(indx),  y0(indy),   z0, v0,1);  

residuals = (obs_rngchg - (modelout));

plot_model(residuals)