function data = GetSpikeWaveforms(unit,varargin)

%GetSpikeWaveforms - Get spike waveforms.
%
%  USAGE
%
%    waveforms = GetSpikeWaveforms(unit,<options>)
%
%    unit           [electrode group,cluster] pair
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'restrict'    list of time intervals to read from the data file
%    =========================================================================
%
%  OUTPUT
%
%    waveforms      3D array (spike #,channel,sample) of waveforms
%
%  SEE
%
%    See also LoadSpikeWaveforms, PlotSpikeWaveforms.

% Copyright (C) 2004-2014 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

% Default values
intervals = [];

if nargin < 1,
	error('Incorrect number of parameters (type ''help <a href="matlab:help GetSpikeWaveforms">GetSpikeWaveforms</a>'' for details).');
end

if ~isivector(unit,'#2','>=0'),
	error('Incorrect unit (type ''help <a href="matlab:help GetSpikeWaveforms">GetSpikeWaveforms</a>'' for details).');
end
group = unit(1);
cluster = unit(2);

% Parse parameter list
for i = 1:2:length(varargin),
  if ~ischar(varargin{i}),
    error(['Parameter ' num2str(i+1) ' is not a property (type ''help <a href="matlab:help GetSpikeWaveforms">GetSpikeWaveforms</a>'' for details).']);
  end
  switch(lower(varargin{i})),
    case {'intervals','restrict'},
      intervals = varargin{i+1};
      if ~isdmatrix(intervals) || size(intervals,2) ~= 2,
			error('Incorrect value for property ''intervals'' (type ''help <a href="matlab:help GetSpikeWaveforms">GetSpikeWaveforms</a>'' for details).');
      end
    otherwise,
		error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help GetSpikeWaveforms">GetSpikeWaveforms</a>'' for details).']);
  end
end

filename = [DATA.session.path '/' DATA.session.basename '.spk.' int2str(group)];
nChannels = length(DATA.spikeGroups.groups{group});
nSamplesPerWaveform = DATA.spikeGroups.nSamples(group);

% Determine which spikes belong to this cluster
t = GetSpikes([group -3],'output','full');
c = t(:,3) == cluster;
if isempty(intervals),
	% List all spikes in this cluster
	list = find(c);
else
	% List spikes that fall in the requested time intervals
	in = InIntervals(t(:,1),intervals);
	list = find(in&c);
end

data = LoadSpikeWaveforms(filename,nChannels,nSamplesPerWaveform,list);
