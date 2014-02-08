%
% sumatra_modes_template.m
% 
% Template script for analyzing modes spectra for sumatra
%
% For a view of the spectra, open this file:
%    /home/admin/databases/SUMATRA/data/wfobject/all_sumatra_modes.pdf
%
% This assumes that you have added the path to the GEOTOOLS directories.
%

clear
close all
clc

iload = 0;      % USER: CHANGE THIS

% add path in order to access additional files and scripts
ddir = '/home/admin/databases/SUMATRA/data/wfobject/';
addpath(ddir);

% first load all the data
if iload==1
    sumatra_modes_fft;
    disp('bad records that will not be used:');
    stas(scut)
    
    % plot the time series that were manually cut from the analysis
    for ii = 1:nw   % loop over all time series
        if any(ii==scut)  
            stag = [stas{ii} '_' chans{ii} '_' nets{ii}];
            stdur = sprintf('duration = %.2f days',get(w(ii),'duration'));
            stit = [stag ', ' stdur];
            disp(sprintf('%i/%i %s',ii,nw,stag));
            figure; hold on; plot(w(ii)); title(stit);
        end
    end
    
    break
end

% load list of files
% note: column 1 is the index into the full set of waveforms [1:169]
%       column 2 is the index into the reduced set of FFT waveforms;
%                this also represents the page number of all_sumatra_modes.pdf
[ind,ind_pdf,sta,chan,net,tag,ikeep] = textread([ddir 'sumatra_modes.txt'],'%f%f%s%s%s%s%f');
for ii=1:length(ind)
   disp(sprintf('%3i %3i %7s %7s %4s %4i',ind(ii),ind_pdf(ii),sta{ii},chan{ii},net{ii},ikeep(ii)));
end
% we only want to consider files with FFTs pre-computed
sta  = sta(logical(ikeep));
chan = chan(logical(ikeep));
net  = net(logical(ikeep));
tag  = tag(logical(ikeep));

% USER: PICK A SET TO PLOT AND SAVE FOR ANALYSIS
% note: these are the same indices as the page numbers of the composite PDF
%ipick = [1:length(sta)];    % view all 139
ipick = [1:3 21];

npick = length(ipick);
w(1,npick) = waveform;  % initialize array of waveforms
pamp = zeros(npick,1);
for ii=1:length(ipick)
    jj = ipick(ii);
    if ikeep(jj)==0
        disp('time series was visibly problematic so spectrum was not computed');
    else
        fname = strcat('w',tag{jj},'.mat');
        ifile = [ddir 'full_length/' fname];
        load(ifile);
        w(ii) = v;

        f = get(v,'fft_freq');
        A = get(v,'fft_amp');
        figure; plot(f*1e3,A); xlim([0.2 1.0]);
        title(tag{jj},'interpreter','none');
        xlabel('frequency, mHz'); ylabel('amplitude');
        
        % save amplitude of a particular peak
        f1 = 0.94; f2 = 0.95;   % CHANGE THESE
        pamp(ii) = max(A(and(f*1e3 > f1, f*1e3 < f2)));
    end
end

% diplay the properties of the object
w(4)
% get properties of the object
[slon,slat,dist_deg,az] = getm(w(4),'STLA','STLO','GCARC','AZ');

% START MODES PROBLEM HERE


%==========================================================================
