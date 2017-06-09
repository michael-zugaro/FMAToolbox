%GlobalSettings - Initialize settings for FMAToolbox.
%
% The global variable SETTINGS holds a number of basic settings such as the
% minimum and maximum distance allowed between LEDs. Modify this variable to
% customize individual values. For instance, to change the maximum distance
% between LEDs, add the following two lines to your startup.m file:
%
%   global SETTINGS;
%   SETTINGS.maxDistance = 100;
%
% To change the global settings, edit <a href="matlab:edit GlobalSettings.m">GlobalSettings.m</a>. Note that this will
% affect all users of the FMAToolbox. User settings take precedence over
% global settings.
%
% For some functions of the toolbox, default values can also be customized
% (see <a href="matlab:help CustomDefaults.m">CustomDefaults.m</a> for details).
%
% SEE
%
%   See also CustomDefaults.

% Copyright (C) 2004-2017 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

global SETTINGS;

if ~isfield(SETTINGS,'minFieldSize'),
	SETTINGS.minFieldSize = 100;
end
