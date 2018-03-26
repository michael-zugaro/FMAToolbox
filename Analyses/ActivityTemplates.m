function [templates,correlations,eigenvalues,eigenvectors] = ActivityTemplates(spikes,varargin)

%ActivityTemplates - Compute activity templates from PCA of spike trains.
%
% Computes the templates for the component activation analysis described in
% Peyrache et al (2009). These templates can then be tested on different data
% sets using <a href="matlab:help ReactivationStrength">ReactivationStrength</a>. Time bins can be automatically determined
% using a fixed bin size and a step, or provided as an explicit list (e.g.
% computed using theta phases).
%
%  USAGE
%
%    [templates,correlations,eigenvalues] = ActivityTemplates(spikes,<options>)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'bins'        list of [start stop] for all bins
%     'binSize'     bin size in s (default = 0.050)
%     'step'        step in s (default = 0.050)
%    =========================================================================
%
%  OUTPUT
%
%    templates      3D array of template matrices (dimension 3 is template #,
%                   ordered in descending order of corresponding eigenvalue)
%    correlations   spike count correlation matrix
%    eigenvalues    significant eigenvalues, listed in descending order
%
%  SEE
%
%    See also ReactivationStrength.

% Copyright (C) 2016-2018 by Michaël Zugaro, Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
bins = [];
defaultBinSize = 0.050;
step = 0.05;
binSize = [];

% Check number of parameters
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
end
% Check parameter sizes
if ~isdmatrix(spikes,'@2'),
	error('Parameter ''spikes'' is not a Nx2 matrix (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
end
% Parse parameter list
for i = 1:2:length(varargin),
    if ~ischar(varargin{i}),
        error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).']);
    end
    switch(lower(varargin{i})),
        case 'binsize',
            binSize = varargin{i+1};
            if ~isdscalar(binSize,'>0'),
                error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
            end
        case 'step',
            step = varargin{i+1};
            if ~isdscalar(step,'>0'),
                error('Incorrect value for property ''step'' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
            end
        case 'bins',
			bins = varargin{i+1};
			if ~isdmatrix(bins,'@2'),
				error('Incorrect value for property ''bins'' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).']);
		end
end

% Options binSize and bins are incompatible
if ~isempty(binSize) && ~isempty(bins),
	error('Parameters ''binSize'' and ''bins'' are incompatible (type ''help <a href="matlab:help ActivityTemplates">ActivityTemplates</a>'' for details).');
end
if isempty(binSize) && isempty(bins),
	binSize = defaultBinSize;
end

templates = NaN;
correlations = NaN;
eigenvalues = NaN;

nUnits = max(spikes(:,2));
if isempty(nUnits), return; end

%% Bin spikes
spikes = sortrows(spikes,1);
id = spikes(:,2);

% Shift spike times to start at 0, and list bins unless explicitly provided
if isempty(bins),
	spikes(:,1) = spikes(:,1) - spikes(1,1);
	bins = (0:step:(spikes(end,1)-binSize))';
	bins(:,2) = bins+binSize;
else
	m = min([min(spikes(:,1)) min(bins(:))]);
	spikes(:,1) = spikes(:,1) - m;
	bins = bins - m;
end

% Create spike count matrix
nBins = size(bins,1);
if isempty(nBins), return; end
n = zeros(nBins,nUnits);
for unit = 1:nUnits,
	n(:,unit) = CountInIntervals(spikes(id==unit,1),bins);
end

%% Create correlation matrix
n = zscore(n);
correlations = (1/(nBins-1))*n'*n;

% Compute eigenvalues/vectors and sort in descending order
[eigenvectors,eigenvalues] = eig(correlations);
[eigenvalues,i] = sort(diag(eigenvalues),'descend');
eigenvectors = eigenvectors(:,i);

%% Keep only significant eigenvalues and compute templates

q = nBins/nUnits;
if q < 1,
	warning('Not enough time bins to determine significant templates');
	eigenvalues = NaN;
end
lambdaMax = (1+sqrt(1/q))^2;
significant = eigenvalues>lambdaMax;
eigenvectors = eigenvectors(:,significant);
templates = zeros(nUnits,nUnits,sum(significant));
for i = 1:sum(significant),
	templates(:,:,i) = eigenvectors(:,i)*eigenvectors(:,i)';
	templates(:,:,i) = templates(:,:,i) - diag(diag(templates(:,:,i))); % remove the diagonal
end
