function plot_covsamples(msamples,rho,tlab,msamples2,rho2,tlab2,mlabs)
% PLOT_COVMSAMPLES generates plots for samples of covariance matrix
%
% INPUT
%    msamples   n x s matrix of vector samples
%    rho        n x n 'analytical' correlation matrix
%    tlab       label for plot
%    msamples2  optional: 2nd set of samples ([] for none)
%    rho2       optional: 2nd 'analytical' correlation matrix ([] for none)
%    tlab2      optional: label for plot
%    mlabs      optional: labels for each variable ([] for default)
%
% EXAMPLE: 
%    plot_covsamples(mpost_samples,rho_post,'mpost',[],[],[],mlabs);
%
% NOTE: We could alternatively estimate the covariance matrix
% (and correlation matrix) directly from the input samples.
%
% called by hw_optim.m
%
% Carl Tape, 01-Jan-2012
%

[n,s] = size(msamples);
disp(sprintf('plot_covsamples.m: n = %i, s = %i',n,s));

if isempty(mlabs)
    %mlabs = strtrim(cellstr(num2str([1:n]')));
    mlabs = repmat(cellstr(''),n,1);
    for ii=1:n, mlabs{ii} = sprintf('i%i',ii); end
end

% whether to plot a second set of samples
if and(~isempty(msamples2),~isempty(rho2)), isecond = 1; else isecond = 0; end

% correlation matrices
% note: we could replace imagesc with a non-toolbox function (pcolor)
figure; nr=1+isecond; nc=1;
subplot(nr,nc,1); imagesc(rho); caxis([-1 1]), colorbar
set(gca,'xtick',[1:n],'xticklabel',mlabs,'xaxislocation','top');
set(gca,'ytick',[1:n],'yticklabel',mlabs);
title(['correlation matrix for ' tlab]); axis equal, axis tight

if isecond==1
subplot(nr,nc,2); imagesc(rho2); caxis([-1 1]), colorbar
set(gca,'xtick',[1:n],'xticklabel',mlabs,'xaxislocation','top');
set(gca,'ytick',[1:n],'yticklabel',mlabs);
title(['correlation matrix for ' tlab2]); axis equal, axis tight
end

% scatterplots
if n > 6
    disp(sprintf('n = %i is > 6, so no scatterplots made',n));
else
    figure; nr=n-1; nc=n-1;
    for ii=1:n-1
       for jj=ii+1:n
           iplot = nc*(ii-1) + jj-1;
           %disp([ii jj iplot]);
           subplot(nr,nc,iplot); hold on;
           plot(msamples(ii,:),msamples(jj,:),'.');
           xlabel(mlabs{ii}); ylabel(mlabs{jj});
           st1 = sprintf('corr(%s) = %.3f',tlab,rho(ii,jj));
           if isecond==1
               plot(msamples2(ii,:),msamples2(jj,:),'r.');
               st2 = sprintf('corr(%s) = %.3f',tlab2,rho2(ii,jj));
               title({st1,st2});
           else
               title(st1);
           end
       end
    end
end

%==========================================================================