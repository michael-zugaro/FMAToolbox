function [channels,groups] = GetChannels(groups)

%GetChannels - Get list of channels for one or more spike groups.
%
%  USAGE
%
%    channels = GetChannels(groups)
%
%    groups         one or more optional electrode groups (default = all groups)
%
%  EXAMPLES
%
%    all = GetChannels;           % list all channels
%    channels = GetChannels(2);   % list all channels for electrode group 2


% Copyright (C) 2004-2018 by Michaël Zugaro, 2018 by Ralitsa Todorova
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global DATA;
if isempty(DATA),
	error('No session defined (did you forget to call SetCurrentSession? Type ''help <a href="matlab:help Data">Data</a>'' for details).');
end

nGroups = DATA.spikeGroups.nGroups;
if nargin < 1,
	groups = 1:nGroups;
else
	if ~isivector(groups,'>0',['<=' int2str(nGroups)]),
		error('Incorrect group list (type ''help <a href="matlab:help GetChannels">GetChannels</a>'' for details).');
	end
end

channels = [DATA.spikeGroups.groups{groups}]';

% Provide information about the group of each channel
if nargout > 1,
    nChannelsPerGroup = cellfun(@length,DATA.spikeGroups.groups)'; 
    nChannelsPerGroup = nChannelsPerGroup(groups);
    if ~isempty(which('repelem')),
        groups = repelem(groups',nChannelsPerGroup);
    else
        g = [];
        for i = 1:length(groups), g = [g; ones(nChannelsPerGroup(i),1)*groups(i)]; end
        groups = g;
    end
end

[channels,order] = sort(channels);
groups = groups(order);