%CCG - Compute multiple cross- and auto-correlograms
%
%  USAGE
%
%    [ccg,t] = CCG(times,id,<options>)
%
%    times          times of all events (sorted)
%    id             ID for each event (e.g. unit ID)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'binSize'     bin size in s (default = 0.01)
%     'duration'    duration in s of each xcorrelogram (default = 2)
%     'smooth'      smoothing size in bins (0 = no smoothing, default)
%     'groups'      group number (1 or 2) for each event, used to restrict
%                   cross-correlograms to pairs across two groups of events
%                   (see EXAMPLE #2 below)
%    =========================================================================
%
%  EXAMPLES
%
%    % Auto- and cross-correlograms between all neurons
%    spikes = GetSpikes('output','numbered');
%    [ccg,t] = CCG(spikes(:,1),spikes(:,2));
%
%    % Only tetrode #1 vs tetrode #2 (e.g. mPFC vs HPC neurons)
%    pfc = GetSpikes([1 -1],'output','numbered');
%    hpc = GetSpikes([2 -1],'output','numbered');
%    m = max(pfc(:,2));
%    [spikes,i] = sortrows([pfc(:,1);hpc(:,1)]);
%    ids = [pfc(:,2);hpc(:,2)+m];
%    ids = ids(i);
%    groups = [ones(size(pfc,1));2*ones(size(hpc,1))];
%    groups = groups(i);
%    [ccg,t] = CCG(s,ids,'groups',groups);
%
%    % Between stimulations and MUA spikes
%    spikes = GetSpikes;
%    stimulatios = GetEvents('Stimulation');
%    d = [spikes(:,1) ones(size(spikes,1)) ; stimulations 2*ones(size(stimulations,1))];
%    d = sortrows(d);
%    [ccg,t] = CCG(d(:,1),d(:,2));
%
%  SEE
%
%    See also ShortTimeCCG.

% Copyright (C) 2012-2013 by Michaël Zugaro, Marie Goutierre
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

function [ccg,t] = CCG(times,id,varargin)

% Default values
duration = 2;
binSize = 0.01;
smooth = 0;
groups = [];

% Check parameters
if nargin < 2,
  error('Incorrect number of parameters (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdvector(times),
	error('Parameter ''times'' is not a real-valued vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(id) && ~isdvector(id),
	error('Parameter ''id'' is not a real-valued scalar or vector (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
if ~isdscalar(id) && length(times) ~= length(id),
	error('Parameters ''times'' and ''id'' have different lengths (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
end
id = id(:);
times = times(:);

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help CCG">CCG</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'binsize',
			binSize = varargin{i+1};
			if ~isdscalar(binSize,'>0'),
				error('Incorrect value for property ''binSize'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'duration',
			duration = varargin{i+1};
			if ~isdscalar(duration,'>0'),
				error('Incorrect value for property ''duration'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'smooth',
			smooth = varargin{i+1};
			if ~isdscalar(smooth,'>=0'),
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
		case 'groups',
			groups = varargin{i+1};
			if ~isdvector(groups) && length(times) ~= length(groups)
				error('Incorrect value for property ''smooth'' (type ''help <a href="matlab:help CCG">CCG</a>'' for details).');
			end
	end
end

% Number of id, number of bins, etc.
if length(id) == 1,
	id = ones(length(times),1);
	nIDs = 1;
else
	nIDs = length(unique(id));
end
halfBins = round(duration/binSize/2);
nBins = 2*halfBins+1;
t = (-halfBins:halfBins)'*binSize;

if length(times) <= 1,
	return
end

% Sort events in time and compute CCGs
[times,i] = sort(times);
id = id(i);
counts = double(CCGHeart(times,uint32(id),binSize,uint32(halfBins)));

% Reshape the results
n = max(id);
counts = reshape(counts,[nBins n n]);
if n < nIDs,
	counts(nBins,nIDs,nIDs) = 0;
end

% Restrict the results to inter-group CCGs if requested
if ~isempty(groups),
	group1 = unique(id(groups == 1));
	group2 = unique(id(groups == 2));
	nGroup1 = length(group1);

	ccg = zeros(nBins,nGroup1,(nIDs-nGroup1));

	for i = 1:nGroup1,
		for j = 1:(nIDs-nGroup1),
			  ccg(:,i,j) = Smooth(flipud(counts(:,group1(i),group2(j))),smooth);
		end
	end
else
	ccg = zeros(nBins,nIDs,nIDs);
	% Reshape
	for g1 = 1:nIDs,
		for g2 = g1:nIDs,
			ccg(:,g1,g2) = Smooth(flipud(counts(:,g1,g2)),smooth);
		end
	end
end
