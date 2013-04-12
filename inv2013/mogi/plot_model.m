function plot_model(D)

% plot an interferogram (both as an unwrapped image and wrapped phase image

SAMPLE = 1100;
LINE = 980;
POSTING = 40.0;
HALF_WAVE = 28.3;

load data_mask

xvec = [1:POSTING:1100*POSTING]./1000;
yvec = [1:POSTING:LINE.*POSTING]./1000;

nanclr = 0.5*[1 1 1];

figure;

fvec = D.*data_mask;

subplot(2,1,1);
imagesc(xvec,yvec,fvec);

% http://stackoverflow.com/questions/8481324/contrasting-color-for-nans-in-imagesc
cm = colormap('jet');
amin = min(D(:));
amax = max(D(:));
n = size(cm,1);
dmap = (amax-amin)/n;
% add nan color to colormap
colormap([nanclr; cm]);
% changing color limits
caxis([amin-dmap amax]);
% place a colorbar
hcb = colorbar;
% change Y limit for colorbar to avoid showing NaN color
ylim(hcb,[amin amax])

%caxis([-30 30])
axis equal; axis tight
set(gca,'FontSize',12); h = xlabel('Easting [km]');set(h,'FontSize',14);
h = ylabel('Northing [km]');set(h,'FontSize',14);
h = title('Line of Sight Motion [mm]');set(h,'FontSize',14);

fvec = wrap(D.*data_mask./10)./5.66.*4*pi;

subplot(2,1,2);
imagesc(xvec,yvec,fvec);
colorbar;
axis equal; axis tight
set(gca,'FontSize',12); h = xlabel('Easting [km]');set(h,'FontSize',14);
h = ylabel('Northing [km]');set(h,'FontSize',14);
h = title('Interferogram phase [rad]');set(h,'FontSize',14);


