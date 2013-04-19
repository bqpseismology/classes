function [xvec,yvec] = plot_model(D)
%PLOT_MODEL plot an interferogram both as an unwrapped image and wrapped phase image

SAMPLE = 1100;
LINE = 980;
POSTING = 40.0;
%HALF_WAVE = 28.3;

xvec = [1:POSTING:SAMPLE*POSTING]/1000;
yvec = [1:POSTING:LINE*POSTING]/1000;
yvec = fliplr(yvec);

F1 = D;
F2 = wrap(D/10)/5.66*4*pi;

iimagesc = 0;           % plot with imagesc (=1) or pcolor (=0)
nanclr = 0.6*[1 1 1];   % NaN color (imagesc=0 only)
itwofig = 0;            % =1 for two figures

if itwofig==0, figure; end

for kk=1:2

    if kk==1
        F = D;
        tlab = 'Line of Sight Motion [mm]';
        clims = [-30 30];
    else
        F = wrap(D/10)/5.66*4*pi;
        tlab = 'Interferogram phase [rad]';
        clims = [-1 1]*max(abs(F(:)));
    end
    
    if itwofig==0, subplot(2,1,kk); else figure; end
    if iimagesc==1
        imagesc(xvec,yvec,F,clims);
        set(gca,'ydir','normal');
    else
        pcolor(xvec,yvec,F); shading flat;
        %set(gca,'ydir','reverse');
        caxis(clims);

        % silly Matlab commands to make sure that the colors that are shown in
        % the figure are actually printed to files
        set(gcf,'Color','white');
        set(gca,'Color',nanclr);
        set(gcf,'InvertHardCopy','off');
    end
    colorbar   
    axis equal; axis tight
    set(gca,'FontSize',12);
    xlabel('Easting [km]','FontSize',14);
    ylabel('Northing [km]','FontSize',14);
    title(tlab,'FontSize',14);

end
