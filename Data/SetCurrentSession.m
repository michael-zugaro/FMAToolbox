function SetCurrentSession(varargin)

%SetCurrentSession - Load all data for a given recording session.
%
% Set current session files and read data from disk. Calling SetCurrentSession
% without parameters will display a file selection dialog.
%
%  USAGE
%
%    SetCurrentSession(filename,<options>)
%
%    filename       optional parameter file name; use 'same' to force reload
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'spikes'      load or skip spike files (default = 'on')
%     'verbose'     display progress messages (default = 'on')
%    =========================================================================
%
%  NOTE
%
%    If no parameter file name is specified, an interactive file selection
%    dialog is displayed.

% Copyright (C) 2004-2017 by Michaël Zugaro, 2014 by Gabrielle Girardeau
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
spikes = 'on';
verbose = 'on';
filename = '';

% Filename?
if nargin ~= 0,
	if ~isastring(varargin{1},'spikes'),
		filename = varargin{1};
		varargin = {varargin{2:end}};
	end
end

% Check number of parameters
if mod(length(varargin),2) ~= 0,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SetCurrentSession">SetCurrentSession</a>'' for details).');
end

% Parse parameter list
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SetCurrentSession">SetCurrentSession</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'verbose',
			verbose = lower(varargin{i+1});
			if ~isastring(verbose,'on','off'),
				error('Incorrect value for property ''verbose'' (type ''help <a href="matlab:help SetCurrentSession">SetCurrentSession</a>'' for details).');
			end
		case 'spikes',
			spikes = lower(varargin{i+1});
			if ~isastring(spikes,'on','off'),
				error('Incorrect value for property ''spikes'' (type ''help <a href="matlab:help SetCurrentSession">SetCurrentSession</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SetCurrentSession">SetCurrentSession</a>'' for details).']);
	end
end

global DATA;
separator = filesep;

% Initialization
if isempty(DATA) || ~isfield(DATA,'session') || ~isfield(DATA.session,'path') || ~isfield(DATA.session,'basename'),
	format long g;
	DATA.session.basename = '';
	DATA.session.path = '';
	DATA.spikeGroups.nGroups = 0;
	DATA.spikeGroups.nSamples = [];
	DATA.spikeGroups.peakSamples = [];
	DATA.spikeGroups.groups = {};
	DATA.nChannels = [];
	DATA.nBits = [];
	DATA.rates.lfp = [];
	DATA.rates.wideband = [];
	DATA.rates.video = [];
	DATA.maxX = [];
	DATA.maxY = [];
	DATA.events.time = [];
	DATA.events.description = {};
	DATA.positions = [];
	DATA.spikes = [];
	% Default settings
	GlobalSettings;
end

if isempty(filename) || (strcmp(filename,'same') && isempty(DATA.session.basename)),
	% Interactive mode
	[filename,path] = uigetfile('*.xml','Please select a parameter file for this session');
	if filename == 0,return; end
	filename = [path filename];
end

if strcmp(filename,'same'),
	% Force reload
	path = DATA.session.path;
	basename = DATA.session.basename;
else
	% Parse file name
	[path,basename] = fileparts(filename);
	if isempty(path),
        path = pwd;
	else
		if ~exist(path),
			error(['Directory ''' path ''' does not exist.']);
		end
		% Clean path (e.g. simplify ../ or ./ substrings) and make it absolute
		here = pwd;
		cd(path);
		path = pwd;
		cd(here);
	end
end

if strcmp(verbose,'on'), disp(['Loading session files for ' basename]); end

% File already loaded?
if strcmp(basename,DATA.session.basename) & strcmp(path,DATA.session.path) & ~strcmp(filename,'same'),
	disp(['... session files already loaded, skipping - type SetCurrentSession(''same'') to force reload']);
	disp('Done');
	return
end

% Parameter file
DATA = LoadParameters([path separator basename '.xml']);
if strcmp(verbose,'on'), disp(['... loaded parameter file ''' basename '.xml''']); end

% Event file(s)
DATA.events.time = [];
DATA.events.description = {};
eventFiles = dir([path separator basename '.*.evt']);
if ~isempty(eventFiles),
	for i = 1:length(eventFiles),
		events = LoadEvents([path separator eventFiles(i).name]);
		if isempty(events.time), continue; end
		DATA.events.time = [DATA.events.time ; events.time];
		DATA.events.description = {DATA.events.description{:} events.description{:}}';
		if strcmp(verbose,'on'), disp(['... loaded event file ''' eventFiles(i).name '''']); end
	end
	[DATA.events.time,i] = sortrows(DATA.events.time);
	DATA.events.description = {DATA.events.description{i}}';
else
	if strcmp(verbose,'on'), disp('... (no event file found)'); end
end

% Position file
DATA.positions = [];
if exist([path separator basename '.pos']),
	DATA.positions = LoadPositions([path separator basename '.pos'],DATA.rates.video);
	if strcmp(verbose,'on'), disp(['... loaded position file ''' basename '.pos''']); end
elseif exist([path separator basename '.whl']),
	DATA.positions = LoadPositions([path separator basename '.whl'],DATA.rates.video);
	if strcmp(verbose,'on'), disp(['... loaded position file ''' basename '.whl''']); end
elseif exist([path separator basename '.whl']),
	DATA.positions = LoadPositions([path separator basename '.mqa'],DATA.rates.video);
	if strcmp(verbose,'on'), disp(['... loaded position file ''' basename '.mqa''']); end
else
	if strcmp(verbose,'on'), disp('... (no position file found)'); end
end

% Spike files
if strcmp(spikes,'on'),
	DATA.spikes = [];
	for i = 1:DATA.spikeGroups.nGroups,
		filename = [path separator basename '.' int2str(i) '.clu'];
		if exist(filename,'file'),
			try
				DATA.spikes = [DATA.spikes;LoadSpikeTimes(filename,DATA.rates.wideband)];
				if strcmp(verbose,'on'), disp(['... loaded spike files ''' basename '.' int2str(i) '.clu''']); end
			catch
				if strcmp(verbose,'on'), disp(['... (could not load spike files ''' basename '.' int2str(i) '.clu'')']); end
			end
		else
			filename = [path separator basename '.clu.' int2str(i)];
			if exist(filename,'file'),
				try
					DATA.spikes = [DATA.spikes;LoadSpikeTimes(filename,DATA.rates.wideband)];
					if strcmp(verbose,'on'), disp(['... loaded spike files ''' basename '.clu.' int2str(i) '''']); end
				catch
					if strcmp(verbose,'on'), disp(['... (could not load spike files ''' basename '.clu.' int2str(i) ''')']); end
				end
			end
		end
	end
	if isempty(DATA.spikes),
		if strcmp(verbose,'on'), disp('... (no spike files found)'); end
	end
else
	if strcmp(verbose,'on'), disp('... (skipping spike files)'); end
end

% This is updated only once the files have been properly loaded
DATA.session.basename = basename;
DATA.session.path = path;

if strcmp(verbose,'on'), disp('Done'); end
