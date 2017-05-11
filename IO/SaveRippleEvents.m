function SaveRippleEvents(filename,ripples,channelID,varargin)

%SaveRippleEvents - Save hippocampal ripple (~200Hz oscillations) events.
%
%  USAGE
%
%    SaveRippleEvents(filename,ripples,channelID,options)
%
%    filename       file to save to
%    ripples        ripple info as provided by <a href="matlab:help FindRipples">FindRipples</a>
%    channelID      channel ID (appended to the event description)
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%     'overwrite'   overwrite file if it exists (default = 'off')
%    =========================================================================
%
%  SEE
%
%    See also FindRipples, RippleStats, PlotRippleStats, SaveEvents.

% Copyright (C) 2004-2015 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% Default values
overwrite = 'off';

if nargin < 3,
  error('Incorrect number of parameters (type ''help <a href="matlab:help SaveRippleEvents">SaveRippleEvents</a>'' for details).');
end

for i = 1:2:length(varargin),
	if ~ischar(varargin{i}),
		error(['Parameter ' num2str(i+2) ' is not a property (type ''help <a href="matlab:help SaveRippleEvents">SaveRippleEvents</a>'' for details).']);
	end
	switch(lower(varargin{i})),
		case 'overwrite',
			overwrite = varargin{i+1};
			if ~isastring(overwrite,'on','off'),
				error('Incorrect value for property ''overwrite'' (type ''help <a href="matlab:help SaveRippleEvents">SaveRippleEvents</a>'' for details).');
			end
		otherwise,
			error(['Unknown property ''' num2str(varargin{i}) ''' (type ''help <a href="matlab:help SaveRippleEvents">SaveRippleEvents</a>'' for details).']);
	end
end

n = size(ripples,1);
r = ripples(:,1:3)';
events.time = r(:);
for i = 1:3:3*n,
	events.description{i,1} = ['Ripple start ' int2str(channelID)];
	events.description{i+1,1} = ['Ripple peak ' int2str(channelID)];
	events.description{i+2,1} = ['Ripple stop ' int2str(channelID)];
end

if strcmp('overwrite','on'),
    SaveEvents(filename,events,'overwrite','on');
else
    SaveEvents(filename,events);
end
