%
% spshell_template.m
% Applied Seismology, Geos 626, UAF
%
% This program computes the toroidal modes for a spherical shell with
% uniform rigidity and density.
%
% calls these functions:
%    stress_disp_tor.m
%    surf_stress.m
%
% GLOBAL VARIABLES: l radius WT rspan mu rho
%
% by Charles Ammon, Penn State, 2000
% modifications by Carl Tape, UAF, 01/2012
%

close all; clear all; clc

% global variables
global l radius WT rspan mu rho % omega

%------------------------------------------------
% USER INPUT

% shell spans from CMB to surface
earthr = 6371000;
rspan = [3480000 earthr];

% Earth model (uniform mantle shell)
% WARNING: BECAUSE THESE ARE GLOBAL VARIABLES, THESE VALUES MAY BE
% OVER-RIDDEN BY VALUES DEFINED IN OTHER FUNCTIONS.
rho = 4380;             % density
mu  = 5930*5930*rho;    % rigidity (mu = 1.54e11 Pa)

l = 2;                  % degree (l >= 1)
rmax = 9;               % maximum number of roots/eigenfunctions/subplots (default = 9)
iplot_eig_freqs = 1;    % plot eigenfunctions (=1) or not (=0)
iplot_all_freqs = 1;

% path to the directory containing the data file prem_Tmodes.txt
ddir = './data/';

iprint = 0; % print figures to file (=1) or not (=0)
pdir = './';

%------------------------------------------------

% range of frequencies (note: omega = 2*pi*f), in Hz
fmin = 1/3600;      % initial frequency to start (T = one hour)
df = 0.0002;        % frequency step size (chosen by trial and error)
%fmax = 0.08;        % stopping frequency (somewhat arbitrary)
fmax = 0.003;
fvec = [fmin:df:fmax];
disp(sprintf('frequency vector ranges from %.3f mHz to %.3f mHz',fmin*1e3,fmax*1e3));
disp(sprintf('num frequency points is %i, df = %.3f mHz',length(fvec),df*1e3));
disp(sprintf('--> period ranges from %.2f min to %.2f min',1/fmin/60,1/fmax/60));

% THIS BLOCK IS FOR PLOTTING ONLY
if iplot_all_freqs==1
    for ii=1:length(fvec)
        % update W(r) and T(r), stored within WT
        surf_stress(fvec(ii));
        
        % plotting parameters
        rplot = radius/1000;
        xmx = 1.1; ymn = rspan(1)/1000; ymx = rspan(2)/1000; dy = 100;

        % displacement for each frequency
        Wplot = WT(:,1)/max(abs(WT(:,1)));
        figure(2); hold on; plot(Wplot,rplot,'b');
        text(Wplot(end),rplot(end)+dy/2,num2str(ii));
        plot([0 0],rspan/1000,'k','linewidth',2);
        xlabel('normalized displacement, W(r)'); ylabel('radius, km');
        axis([-xmx xmx ymn-dy ymx+dy]);

        % stress for each frequency
        Tplot = WT(:,2)/max(abs(WT(:,2)));
        figure(3); hold on; plot(Tplot,rplot,'r');
        text(Tplot(end),rplot(end)+dy/2,num2str(ii));
        plot([0 0],rspan/1000,'k','linewidth',2);
        xlabel('normalized stress, T(r)'); ylabel('radius, km');
        axis([-xmx xmx ymn-dy ymx+dy]);
    end
    
    % print figures for HW
    if iprint==1
        figure(2); print(gcf,'-depsc',[pdir 'modes_Wr']);
        figure(3); print(gcf,'-depsc',[pdir 'modes_Tr']);
    end

    % exit
    break
end

% initial freqeuncy and corresponding surface stress
f = fvec(1);
Tsurf = surf_stress(f);
nn = 1;             % counter for the root number (starting at one)

% leave gap for T(n=0,l=1), which do not exist
% note: these are useful when looping over degree l
if and(l==1,rmax>1), nn = 2; end        % fill the n >= 1 entries
if and(l==1,rmax==1), continue; end     % exit loop early

% THIS IS THE KEY LOOP OVER FREQUENCIES
froots = NaN*ones(rmax,1);
for ii = 2:length(fvec)-1
    % frequency interval over which we check for a root
    oldf = f;
    f = fvec(ii);
    
    % The function surf_stress.m will updated the key variable WT,
    % which contains the radial displacement W(r) in the first column
    % and stress T(r) in the second column.
    Tsurfold = Tsurf;          % surface stress for previous f
    Tsurf = surf_stress(f);    % surface stress for new f
    
    disp(sprintf('%3i %10.3e %10.3e %.2f mHz %.1f s %.2f min', ...
        ii, Tsurfold, Tsurf, f*1e3, 1/f, 1/f/60));

    % Check if the value of the surface-stress has changed sign,
    % which would indicate that we passed at least one root.
    % If we did cross a root, call the matlab function fzero to refine the root.
    % Then store the root in the vector froots and plot the results.
    if (Tsurfold * Tsurf < 0)
        f0 = fzero('surf_stress',[oldf f]);
        froots(nn) = f0;

        % update eigenfunctions (WT, radius) for the exact frequency
        surf_stress(f0);
        disp(sprintf('%3i %3i %.2f mHz', ii, nn-1, f0*1e3));

        % plotting eigenfunctions (displacement and stress as a function of radius)
        if iplot_eig_freqs==1
            xmx = 1.2; ymn = rspan(1)/1000; ymx = rspan(2)/1000;
            figure(1); if rmax==1, subplot(1,1,nn); else subplot(3,3,nn); end
            hold on;
            plot(WT(:,1)/max(abs(WT(:,1))),radius/1000,'b');    % W(r), displacement (blue)
            plot(WT(:,2)/max(abs(WT(:,2))),radius/1000,'r');    % T(r), stress (red)
            plot([-xmx xmx],ymn*[1 1],'k',[-xmx xmx],ymx*[1 1],'k',[0 0],[ymn ymx],'k');
            axis([-xmx xmx ymn-300 ymx+300]); %grid on;
            title(sprintf('f = %.2f mHz, T = %.2f min (l = %i)',froots(nn)*1000,1/f0/60,l));
            text(-1,ymx+100,sprintf('n = %i',nn-1));
            if mod(nn-1,3)==0, ylabel('radius (km)'); end
        end
        
        if nn==rmax, break; end % this halts the search after rmax roots
        nn = nn + 1;
    end
end
fprintf('l = %i, nroot = %i (rmax = %i)\n',l,sum(~isnan(froots)),rmax);

break

% observations used in PREM
dfile = [ddir 'prem_Tmodes.txt'];
[nobs,~,lobs,T,Tstd] = textread(dfile,'%f%s%f%f%f','headerlines',6);
disp('normal mode observations (measured from seismograms):');
for ii=1:length(nobs)
   disp(sprintf('n = %i, l = %2i, T = %8.2f +/- %.2f s',nobs(ii),lobs(ii),T(ii),Tstd(ii))); 
end

%==========================================================================
