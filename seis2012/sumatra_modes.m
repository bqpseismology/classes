% sumatra_modes.m
% 
% Template script for analyzing modes spectra for sumatra
%
% For a view of the spectra, open this file:
%    /home/admin/databases/SUMATRA/data/wfobject/all_sumatra_modes.pdf
%

clear
close all
clc

% CHANGE THIS TO YOUR BASE GEOTOOLS DIRECTORY
gdir = '/home/carltape/GEOTOOLS/';

iload = 1;      % CHANGE THIS
ddir = '/home/admin/databases/SUMATRA/data/wfobject/';

% add paths (in principle, these only need to be executed once)
addpath([gdir 'matlab_util/']); % GEOTOOLS (getwaveform.m)
addpath(ddir);                  % sumatra_modes_fft.m

% first load all the data
if iload == 1
    sumatra_modes_fft;

    disp('bad records that will not be used:');
    stas(scut)
    
    % plot the time series that were manually cut from the analysis
    for ii = 1:nw
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
[ind,sta,chan,net,tag,ikeep] = textread([ddir 'sumatra_modes.txt'],'%f%s%s%s%s%f');
for ii=1:length(ind)
   disp(sprintf('%3i %7s %7s %4s %4i',ii,sta{ii},chan{ii},net{ii},ikeep(ii)));
end

% USER: PICK A SET TO PLOT AND SAVE FOR ANALYSIS
ipick = [1:3 22];

npick = length(ipick);
w(1,npick) = waveform;  % initialize array of waveforms
for ii=1:length(ipick)
    jj = ipick(ii);
    if ikeep(jj)==0
        disp('time series was visibly problematic so spectrum was not computed');
    else
        fname = strcat('w',tag{jj},'.mat');
        ifile = [ddir 'full_length/' fname];
        load(fname);
        w(ii) = v;

        f = get(v,'fft_freq');
        A = get(v,'fft_amp');
        figure; plot(f*1e3,A); xlim([0.2 1.0]);
        title(tag{jj},'interpreter','none');
        xlabel('frequency, mHz'); ylabel('amplitude');
    end
end

% diplay the properties of the object
w(4)
[slon,slat,dist_deg,az] = getm(w(4),'STLA','STLO','GCARC','AZ');

% START MODES PROBLEM HERE


%==========================================================================