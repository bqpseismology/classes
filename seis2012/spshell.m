%
% spshell.m
% Applied Seismology, UAF
%
% This program computes the toroidal modes for a spherical shell with
% uniform rigidity and density.
%
% calls these functions:
%    stress_disp_tor.m
%    surf_stress.m
%    earthfun.m
%
% GLOBAL VARIABLES: l radius sd rspan mu rho
%
% by Charles Ammon, Penn State, 2000
% modifications by Carl Tape, UAF, 01/2012
%

close all
clear all
clc

% global variables
global l radius sd rspan mu rho % omega

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

l = 2;         % degree (l >= 1)
rmax = 9;       % maximum number of roots/eigenfunctions/subplots (default = 9)
iploteig = 1;   % plot eigenfunctions (=1) or not (=0)

%------------------------------------------------

% range of frequencies (note: omega = 2*pi*f)
fmin = 1/3600;      % initial frequency to start (T = one hour)
df = 0.0002;        % frequency step size (chosen by trial and error)
fmax = 0.08;        % stopping frequency (arbitrary)
fvec = [fmin:df:fmax];

% initial freqeuncy and corresponding surface stress
f = fvec(1);
surface_stress = surf_stress(f);
jj = 1;             % counter for the root number (starting at one)

% leave gap for T(n=0,l=1), which do not exist
if and(l==1,rmax>1), jj = 2; end

froots = NaN*ones(rmax,1);
for ii = 2:length(fvec)-1   % loop over frequencies
    oldf = f;
    f = fvec(ii);

    % compute the surface stress for this frequency
    % this is calling the scripts listed above
    oldvalue = surface_stress;
    surface_stress = surf_stress(f);

    disp(sprintf('%3i %10.3e %10.3e %.2f mHz %.1f s %.2f min', ...
        ii, oldvalue, surface_stress, f*1e3, 1/f, 1/f/60));

    % Check if the value of the surface-stress has changed sign
    % which would indicate that we passed at least one root.
    % If we did cross a root, call the matlab function fzero to refine the root.
    % Then store the root in the vector root() and plot the results.
    if (oldvalue * surface_stress < 0)
        f0 = fzero('surf_stress',[oldf f]);
        froots(jj) = f0;

        % update eigenfunctions (sd, radius) for the exact frequency
        surf_stress(f0);
        disp(sprintf('%3i %.2f mHz', ii, f0*1e3));

        % plotting eigenfunctions (displacement and stress as a function of radius)
        if iploteig==1
            xmx = 1.2; ymn = rspan(1)/1000; ymx = rspan(2)/1000;
            figure(1); if rmax==1, subplot(1,1,jj); else subplot(3,3,jj); end
            hold on;
            plot(sd(:,1)/max(abs(sd(:,1))),radius/1000,'b');    % displacement (blue)
            plot(sd(:,2)/max(abs(sd(:,2))),radius/1000,'r');    % stress (red)
            plot([-xmx xmx],ymn*[1 1],'k',[-xmx xmx],ymx*[1 1],'k',[0 0],[ymn ymx],'k');
            axis([-xmx xmx ymn-300 ymx+300]); %grid on;
            title(sprintf('f = %.2f mHz, T = %.2f min (l = %i)',froots(jj)*1000,1/f0/60,l));
            text(-1,ymx+100,sprintf('n = %i',jj-1));
            if mod(jj-1,3)==0, ylabel('radius (km)'); end
        end
        
        if jj==rmax, break; end % this halts the search after rmax roots
        jj = jj + 1;
    end
end
fprintf('l = %i, nroot = %i (rmax = %i)\n',l,sum(~isnan(froots)),rmax);

% observations used in PREM (assumes lvec = 1:10)
dfile = './data/prem_Tmodes.txt';
[nobs,~,lobs,T,Tstd] = textread(dfile,'%f%s%f%f%f','headerlines',6);
disp('normal mode observations:');
for ii=1:length(nobs)
   disp(sprintf('n = %i, l = %2i, T = %8.2f +/- %.2f s',nobs(ii),lobs(ii),T(ii),Tstd(ii))); 
end

%==========================================================================
