function output = TimeToPhase(input,phases,varargin)

%TimeToPhase - Resample signal as a function of unwrapped phase.
%
%  LFP signals and spikes are recorded as a function of computer time, but one
%  may instead wish to express them as a function of unwrapped LFP phase. This
%  is necessary e.g. to align LFPs and spikes with reconstructed positions
%  (see <a href="matlab:help ReconstructPosition">ReconstructPosition</a>), computing spike autocorrelograms in units of
%  'per-cycle' rather than Hz, estimating cross-frequency (amplitude-phase)
%  coupling, etc.
%
%  To better understand how this transformation acts on the signal, consider the
%  following example. When plotting LFP traces as a function of time, unless the
%  underlying rhythm (e.g. theta) is a pure sine wave, phases will increase in a
%  non-linear way along the x axis, sometimes faster, sometimes slower. Resampling
%  in phase space 'distorts' the signal and yields a different plot where phases
%  increase linearly (but computer clock implicitly increases non-linearly).
%
%
%  USAGE
%
%    output = TimeToPhase(input,phases,varargin)
%
%    input          input <a href="matlab:help samples">samples</a>
%    phases         unwrapped phases (obtained e.g. using <a href="matlab:help Phase">Phase</a>)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'frequency'   number of output samples per cycle (default = 150)
%    =========================================================================
%

% Copyright (C) 2015 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Defaults
frequency = 150;

% Check number of parameters
if nargin < 2 || mod(length(varargin),2) ~= 0,
	error('Incorrect number of parameters (type ''help <a href="matlab:help TimeToPhase">TimeToPhase</a>'' for details).');
end

% Check parameters
if ~issamples(input),
	error('Incorrect input samples (type ''help <a href="matlab:help TimeToPhase">TimeToPhase</a>'' for details).');
end

% Parse parameters
for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help TimeToPhase">TimeToPhase</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'frequency',
			frequency = varargin{i+1};
			if ~isdscalar(frequency),
				error('Incorrect value for property ''frequency'' (type ''help <a href="matlab:help TimeToPhase">TimeToPhase</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help TimeToPhase">TimeToPhase</a>'' for details).']);
	end
end

if size(input,2) == 1,
	% Point process
	[~,~,output] = Phase(phase,input);
else
	% Continuous data
	start = phases(1,2);
	stop = phases(end,2);
	nCycles = (stop-start)/(2*pi);
	output = Interpolate([phases(:,2) input(:,2)],linspace(start,stop,round(nCycles*frequency)));
end

