function [events,info] = FindSpikeTrains(spikes,varargin);

%FindSpikeTrains - Split long spike trains into 'active' periods.
%
% Split spike trains into 'active' events in which consecutive spikes remain
% in close temporal proximity. These events are separated by periods of
% silence, cannot exceed a certain duration, and must recruit a minimum number
% of units.
%
%  USAGE
%
%    [info,events] = FindSpikeTrains(spikes,<options>)
%
%    spikes         list of (t,id) pairs, e.g. obtained using <a href="matlab:help GetSpikes">GetSpikes</a> with
%                   option 'output' set to 'numbered'
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'maxDuration' maximum event duration (default = 1 s)
%     'gap'         maximum inter-spike interval in an event (default = 0.1 s)
%     'minUnits'    minimum number of units involved in an event (default = 3)
%    =========================================================================
%
%  OUTPUT
%
%    events         list of (start,end) pairs
%    info           list of (event #,t,unit #) triplets

% Copyright (C) 2016 by Ralitsa Todorova, Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
maxDuration = 1;
gap = 0.1;
minUnits = 3;

% Check parameters
if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help FindSpikeTrains">FindSpikeTrains</a>'' for details).');
end
if ~isdmatrix(spikes,'@2'),
  error('Incorrect spikes (type ''help <a href="matlab:help FindSpikeTrains">FindSpikeTrains</a>'' for details).');
end

% Parse options
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help FindSpikeTrains">FindSpikeTrains</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'maxduration',
			maxDuration = varargin{i+1};
			if ~isdscalar(maxDuration,'>0'),
				error('Incorrect value for property ''maxDuration'' (type ''help <a href="matlab:help FindSpikeTrains">FindSpikeTrains</a>'' for details).');
			end
		case 'gap',
			gap = varargin{i+1};
			if ~isdscalar(gap,'>0'),
				error('Incorrect value for property ''gap'' (type ''help <a href="matlab:help FindSpikeTrains">FindSpikeTrains</a>'' for details).');
			end
		case 'minunits',
			minUnits = varargin{i+1};
			if ~isdscalar(minUnits,'>0'),
				error('Incorrect value for property ''minUnits'' (type ''help <a href="matlab:help FindSpikeTrains">FindSpikeTrains</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help FindSpikeTrains">FindSpikeTrains</a>'' for details).']);
	end
end

spikes = sortrows(spikes);

% Find gaps
gaps = [1;(diff(spikes(:,1))>gap)];
starts = strfind(gaps',[1 0])';
ends = strfind(gaps',[0 1])';

% Make sure the last spike train has an end
if length(ends) < length(starts),
	ends(end+1) = length(spikes);
end

% Count units in each event
for i=1:length(starts),
    nUnits(i,1) = numel(unique(spikes(starts(i):ends(i),2)));
end

% Change starts and ends from indices to timestamps
starts = spikes(starts,1);
ends = spikes(ends,1);

% Event selection criteria: n units, duration
enoughUnits = nUnits >= minUnits;
tooLong = ends-starts > maxDuration;

% Keep selected events
events = [starts(enoughUnits&~tooLong) ends(enoughUnits&~tooLong)];

% Determine spike information
[in,which] = InIntervals(spikes,events);
info = [which(in) spikes(in,:)];
