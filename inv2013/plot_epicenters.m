function plot_epicenters(mprior_samples,mprior,minitial,mtarget,opts,mpost)
% PLOT_EPICENTERS plots epicenters within the source-station geometry
%
% called by forward_epicenter.m, optimization.m
%
% Carl Tape, 24-Feb-2010
%

% input options
xrec = opts{1};
yrec = opts{2};
iray = opts{3};
ndata = length(xrec);

figure; hold on;
msizer = 16;    % receiver size
msizes = 10;    % source size
rthick = 1;     % receiver edge thickness
sthick = 2;     % source edge thickness
rfsize = 10;

% plot ray paths
if iray==1
    for ii=1:ndata
        plot([minitial(2) xrec(ii)],[minitial(3) yrec(ii)],'k','linewidth',1);
    end
end

% plot option depends on if a posterior model is passed
if exist('mpost','var')
    p0 = plot(mprior_samples(2,:),mprior_samples(3,:),'.');
    p1 = plot(minitial(2),minitial(3),'o','markersize',msizes,'markerfacecolor','k','markeredgecolor','w','linewidth',sthick);
    p2 = plot(mpost(2),mpost(3),'o','markersize',msizes,'markerfacecolor','c','markeredgecolor','w','linewidth',sthick);
    pP = plot(mprior(2),mprior(3),'o','markersize',msizes,'markerfacecolor','b','markeredgecolor','w','linewidth',sthick);
    pT = plot(mtarget(2),mtarget(3),'o','markersize',msizes,'markerfacecolor','r','markeredgecolor','w','linewidth',sthick);
    plot(xrec,yrec,'kV','markersize',msizer,'linewidth',rthick);
    for ii=1:ndata, text(xrec(ii),yrec(ii),num2str(ii),'fontsize',rfsize,'color','r',...
            'horizontalalignment','center','verticalalignment','middle'); end
    legend([p0(1) p1 p2 pP pT],'Cprior sample','minitial','mpost','mprior','mtarget');
    %legend([p0(1) p1 p2 pP pT],'Cprior sample','m00',sprintf('m%2.2i',niter),'mprior','mtarget');
    
else
    if ~isempty(mprior_samples)
        p0 = plot(mprior_samples(2,:),mprior_samples(3,:),'.');
    end
    p1 = plot(minitial(2),minitial(3),'o','markersize',msizes,'markerfacecolor','k','markeredgecolor','w','linewidth',sthick);
    pP = plot(mprior(2),mprior(3),'o','markersize',msizes,'markerfacecolor','b','markeredgecolor','w','linewidth',sthick);
    pT = plot(mtarget(2),mtarget(3),'o','markersize',msizes,'markerfacecolor','r','markeredgecolor','w','linewidth',sthick);
    plot(xrec,yrec,'kV','markersize',msizer,'linewidth',rthick);
    for ii=1:ndata, text(xrec(ii),yrec(ii),num2str(ii),'fontsize',rfsize,'color','r',...
            'horizontalalignment','center','verticalalignment','middle'); end
    if ~isempty(mprior_samples)
        legend([p0(1) p1 pP pT],'Cprior sample','minitial','mprior','mtarget');
    else
        legend([p1 pP pT],'minitial','mprior','mtarget');
    end
end

axis equal; axis([0 100 0 100]); %grid on;
set(gca,'xtick',[0:20:100],'ytick',[0:20:100]);
xlabel(' X distance (km)'); ylabel(' Y distance (km)');

%==========================================================================