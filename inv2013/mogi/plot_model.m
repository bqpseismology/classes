function plot_model(modeloutput)

% plot an interferogram (both as an unwrapped image and wrapped phase image

SAMPLE = 1100;
LINE = 980;
POSTING = 40.0;
HALF_WAVE = 28.3;
load data_mask
figure;
subplot(2,1,1);imagesc([1:POSTING:1100*POSTING]./1000,[1:POSTING:LINE.*POSTING]./1000,modeloutput.*data_mask);colorbar;caxis([-30 30])
axis equal; axis tight
set(gca,'FontSize',12); h = xlabel('Easting [km]');set(h,'FontSize',14);
h = ylabel('Northing [km]');set(h,'FontSize',14);
h = title('Line of Sight Motion [mm]');set(h,'FontSize',14);

subplot(2,1,2);imagesc([1:POSTING:1100*POSTING]./1000,[1:POSTING:LINE.*POSTING]./1000,wrap(modeloutput.*data_mask./10)./5.66.*4*pi);colorbar;axis equal; axis tight
set(gca,'FontSize',12); h = xlabel('Easting [km]');set(h,'FontSize',14);
h = ylabel('Northing [km]');set(h,'FontSize',14);
h = title('Interferogram phase [rad]');set(h,'FontSize',14);