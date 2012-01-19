function [Ncum,N,Medges] = seis2GR(mag,dmag)
% SEIS2GR convert catalog to binned magnitudes for analysis
%
%   Converts seismicity catalog to Gutenberg-Richter frequency-magnitude
%   distribution, in both a cumulative version and not cumulative version.
%
% Carl Tape, 20-Oct-2011
%
  
% USER INPUT
idisplay = 1;
ifigure = 0;

n = length(mag);
minm = min(mag);
maxm = max(mag);

disp(sprintf('seis2GR.m: %i events, min M = %.3f, max M = %.3f',n,minm,maxm));

emin = floor(minm/dmag)*dmag;
emax = ceil(maxm/dmag)*dmag;
Medges = [emin:dmag:emax]';

% The last bin will count any values of X that match EDGES(end).
% You can typically ignore the final entry of N and Ncum, which are both zero.
[N,bin] = histc(mag,Medges);
Ncum = flipud(cumsum(flipud(N)));       % from right to left
nbin = length(N);

if idisplay == 1
    for ii=1:nbin-1
        disp(sprintf('Mw = [%5.2f, %5.2f], N = %6.0f, cNL = %7.0f',...
            Medges(ii),Medges(ii+1),N(ii),Ncum(ii)));
    end
end

if ifigure == 1
    figure; nr=2; nc=2;
    for kk=1:4
        if kk==1, D = N; ylab = [' Number (N = ' num2str(n) ')']; end
        if kk==2, D = log10(N); ylab = [' log10[Number] (N = ' num2str(n) ')']; end
        if kk==3, D = Ncum; ylab = [' Cumulative Number (N = ' num2str(n) ')']; end
        if kk==4, D = log10(Ncum); ylab = [' log10[Cumulative Number] (N = ' num2str(n) ')']; end
        subplot(nr,nc,kk); bar(Medges,D,'histc');
        xlim([min(Medges) max(Medges)]); xlabel('Magnitude'); ylabel(ylab);
        h = findobj(gca,'Type','patch'); set(h,'FaceColor',[0 1 1],'EdgeColor','k');
    end
end

%==========================================================================
% EXAMPLES

if 0==1
    % AEIC catalog (from UAF Linux network only)
    catalog = aeic_catalog.get_total;
    mag = catalog.prefMagnitude;
    
    DX = 0.2;
    [Ncum,N,Medges] = seis2GR(mag,DX);
    Mcenter = Medges(1:end-1) + DX;
    Ncumcenter = Ncum(1:end-1);
end

%==========================================================================
